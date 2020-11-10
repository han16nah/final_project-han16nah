#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Course      : Advanced Geoscripting 2020
# Date        : 2020/11/10
# Name        : Hannah Weiser
# e-mail      : h.weiser@stud.uni-heidelberg.de
# =============================================================================
"""
This script performs an ecological data analysis of species occurrence. The steps are the following:
- retrieving data from the Eurostat GISCO API and the GBIF API and converting it to GeoDataFrames
- plotting point observations coloured by the year of observations
- counting the number of point observations in administrative units and plotting choropleth maps
- deriving a co-occurrence matrix and plotting the co-occurrence of one species (e.g. a predator) to several other
species (e.g. preys)

usage: python gbif_analysis.py

or uncomment line 181 and do:
python gbif_analysis.py <config.json>
"""
# =============================================================================
# Imports
# =============================================================================
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from pygbif import occurrences as occ  # python client for the GBIF API
from pygbif import species
import sys
import os
import requests
import json
import geojson
plt.style.use('seaborn')  # figure style


def get_nuts(spatialtype='RG', resolution='10M', year='2021', projection='4326', subset='LEVL_3', output_file="NUTS"):
    """
    Downloads NUTS data (statistical administrative data) of europe.
    For parameters, see also: https://gisco-services.ec.europa.eu/distribution/v2/nuts/nuts-2021-files.html
    :param spatialtype: spatial type; default: 'RG'
        'BN' = boundaries (multilines),
        'RG' = regions (multipolygons),
        'LB' = labels (points)
    :param resolution: resolution; default: '10M'
        '60M' = 1 : 60 million
        '20M' = 1 : 20 million
        '10M' = 1 : 10 million
        '03M' = 1 : 3 million
        '01M' = 1 : 1 million
    :param year: the year of NUTS regulation; default: '2021'
        '2021' / '2016' / '2013' / '2010' / '2006' / '2003'
    :param projection: 4-digit EPSG-code; default: '4326'
        '4326' = GS84, coordinates in decimal degrees
        '3035' = ETRS 1989 in Lambert Azimutal projection
        '3857' = WGS84 Web Mercator Auxiliary Sphere, coordinates in meters
    :param subset: default: 'LEVL_3',
        'LEVL_0': NUTS level 0 (countries)
        'LEVL_1': NUTS level 1
        'LEVL_2': NUTS level 2
        'LEVL_3': NUTS level 3
         No subset means all NUTS levels are in the same file
    :param output_file: Relative or absolute path of output file
    :return: GeoDataFrame of downloaded NUTS data
    """
    theme = 'NUTS'
    file_format = 'geojson'
    file = "%s_%s_%s_%s_%s_%s.%s" % (theme, spatialtype, resolution, year, projection, subset, file_format)
    response = requests.get("https://gisco-services.ec.europa.eu/distribution/v2/nuts/geojson/" + file)
    filename, ext = os.path.splitext(output_file)
    if ext == "":
        output_file += ".geojson"
    elif ext != ".geojson":
        output_file = filename + ".geojson"
        print("Bad file format provided. Changed file extension to .geojson")
    with open(output_file, 'w') as f:
        geojson.dump(response.json(), f)
    gdf = gpd.read_file(output_file)

    return gdf


def occs_to_gdf(occurrence_dict):
    """
    Converts an occurrence dictionary to a geopandas GeoDataFrame and projects it to WGS84 (crs of occurrences in GBIF).
    :param occurrence_dict: occurrence dictionary as returned by the GBIF API (occurrence section)
    :return: occurrence GeoDataFrame
    """
    occ_df = pd.DataFrame.from_dict(occurrence_dict['results'])
    occ_df = occ_df.dropna(subset=['decimalLatitude', 'decimalLongitude'])
    occ_df = occ_df.astype({'year': 'int'})
    occurrence_gdf = gpd.GeoDataFrame(occ_df,
                                      geometry=gpd.points_from_xy(occ_df.decimalLongitude, occ_df.decimalLatitude))
    occurrence_gdf.set_crs('EPSG:4326', inplace=True)

    return occurrence_gdf


def clip_gdf(gdf, clip_feat, crs=None):
    """
    Clips a DataFrame by a clip feature and makes sure they are both in the same crs.
    :param gdf: Input GeoDataFrame to clip.
    :param clip_feat: Clip feature (GeoDataFrame or GeoSeries)
    :param crs: Coordinate reference system to project the input data to, in the form "EPSG:4326"
        If "None", will use the crs of gdf
    :return: Clipped GeoDataFrame
    """
    # Catch the error if one of the input DataFrames has no coordinate reference system
    try:
        if crs is None:
            crs = gdf.crs
            clip_feat.to_crs(crs, inplace=True)
        else:
            clip_feat.to_crs(crs, inplace=True)
            gdf.to_crs(crs, inplace=True)
    except Exception as err:
        print("ValueError: GeoDataFrame has no coordinate reference system. Please set it first.")
        print("Geopandas err message:", err)
        sys.exit()
    gdf_clipped = gpd.clip(gdf, clip_feat)

    return gdf_clipped


def sjoin_and_merge(polys, points, pivot_index, count_column):
    """
    Performs a spatial join and a merge to count the number of "points" with specific attributes
     falling into each polygon item in "polys"
    :param polys: GeoDataFrame containing (Multi-)Polygon geometries
    :param points: GeoDataFrame containing points
    :param pivot_index: Index for creation of the pivot table
    :param count_column: Column values to count
    :return: poly GeoDataFrame with new columns, named by values of count_column and containing the column counts
    """
    try:
        crs = polys.crs
    except Exception as err:
        print("ValueError: 'Polys' GeoDataFrame has no coordinate reference system. Please set it first.")
        print("Geopandas err message:", err)
        sys.exit()
    try:
        points.set_crs(crs, inplace=True)
    except Exception as err:
        print("ValueError: 'Points' GeoDataFrame has no coordinate reference system. Please set it first.")
        print("Geopandas err message:", err)
        sys.exit()
    dfsjoin = gpd.sjoin(polys, points)
    dfpivot = pd.pivot_table(dfsjoin, index=pivot_index, columns=count_column, aggfunc={count_column: len})
    dfpivot.columns = dfpivot.columns.droplevel()
    polys_new = polys.merge(dfpivot, how='left', on=pivot_index)

    return polys_new


def build_cooccurrence_matrix(df):
    """
    Reclassifies the values in a dataframe and then builds  a co-occurrence matrix.
    Values > 0      are reclassified to     1
    NaN-values      are reclassified to     0
    :param df: DataFrame to build the co-occurrence matrix for
    :return: Tuple:
        (Co-occurrence matrix, Co-occurrence matrix with diagonals filled with 0)
    """
    # reclassify: values > 0 --> 1; NaN --> 0
    df = df.where(~(df.values > 0), other=1)
    df = df.fillna(0)
    df_asint = df.astype(int)
    df_cooc = df_asint.T.dot(df_asint)
    df_cooc_f = df_cooc.copy()
    np.fill_diagonal(df_cooc_f.values, 0)

    return df_cooc, df_cooc_f


print("Getting inputs from config file ...")
# input data is read from config-file
config_file = "config.json"
# uncomment following line to provide configuration file in command line
# config_file = sys.argv[1]
with open(config_file) as src:
    config = json.load(src)
try:
    species_list = config["species_list"]
except KeyError:
    print("No species list provided. Aborting.")
    sys.exit()
try:
    country_code = config["country_code"]
except KeyError:
    country_code = None
    print("Key not found. 'country_code' was set to None.")
try:
    out_dir = config["output_dir"]
except KeyError:
    out_dir = './output'
    print("Key not found. 'output_dir' was set to './output'.")
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

print("Downloading NUTS data ...")
# Download NUTS data from Eurostat GISCO API and save to gdf; we use EPSG=3035 because we want an equal-area projection
gdf_aoi_adm = get_nuts(projection="3035")
# if a country_code is provided, extract the subset from only this country
if country_code is not None:
    gdf_aoi_adm = gdf_aoi_adm[gdf_aoi_adm['CNTR_CODE'] == country_code]
gdf_aoi_occs = gdf_aoi_adm.copy()

print("Getting species keys ...")
# get the species keys for the provided species using the species module of the GBIF API
species_keys = [species.name_backbone(x)['usageKey'] for x in species_list]

# iterate over each species, download the data from the GBIF API, write to csv, convert to GeoDataFrame and plot
# the occurrence coloured by the year of observation.
# Finally, add the observation count of each species in each administrative unit to the NUTS GeoDataFrame
# in new columns named after the taxon key
for i, sp in enumerate(species_list):
    print("Processing species nr. %s ..." % (i+1))
    occ_dict = occ.search(scientificName=sp, country=country_code)
    if not occ_dict['results']:
        continue
    occ_gdf = occs_to_gdf(occ_dict)
    # project to the same crs as exported NUTS data ('EPSG:3035')
    occ_gdf.to_crs(epsg=3035, inplace=True)

    occ_gdf.to_csv(os.path.join(out_dir, "gbif_occ_" + occ_dict['results'][0]['scientificName'].
                                replace(" ", "_") + ".csv"))
    fig, ax = plt.subplots(figsize=(10, 10))
    gdf_aoi_adm.plot(ax=ax, color='lightgrey')
    gdf_aoi_adm.geometry.boundary.plot(ax=ax, alpha=0.25, color='black', lw=0.5)
    occ_gdf.plot(column="year", ax=ax, cmap="OrRd", legend=True, categorical=True)
    plt.title("%s - Observations by year" % sp)
    plt.savefig(os.path.join(out_dir, sp + "_occ_by_year.png"))
    plt.clf()
    print("Saved figure to " + os.path.join(out_dir, sp + "_occ_by_year.png"))
    gdf_aoi_occs = sjoin_and_merge(gdf_aoi_occs, occ_gdf, 'NUTS_ID', 'taxonKey')

print("Plotting choropleth map ...")
# Plot choropleth maps showing the number of observations in each administrative unit
n_plots = len([key for key in species_keys if key in gdf_aoi_occs.columns])
n_cols = 4
n_rows = int(np.ceil(n_plots/n_cols))
fig, axs = plt.subplots(n_rows, n_cols, figsize=(20, 15))
i = 0
for sp in species_list:
    col = species_keys[species_list.index(sp)]
    if col not in gdf_aoi_occs.columns:
        continue
    plt_row = int(np.floor(i/n_cols))
    plt_col = i-plt_row*n_cols
    ax = axs[plt_row][plt_col]
    gdf_aoi_occs.plot(ax=ax, color="lightgrey")
    gdf_aoi_occs.plot(column=col, cmap='OrRd', ax=ax, legend=True)
    ax.title.set_text(sp)
    i += 1
plt.savefig(os.path.join(out_dir, "choropleth.png"))
print("Saved figure to " + os.path.join(out_dir, "choropleth.png"))
plt.clf()

# Build a co-occurrence matrix of all species
keys = ['NUTS_ID'] + species_keys
keys = [key for key in keys if key in gdf_aoi_occs.columns]
df_aoi_occs_sub = gdf_aoi_occs[keys]
df_aoi_occs_sub.set_index('NUTS_ID', inplace=True)
df_aoi_cooc, df_aoi_cooc_f = build_cooccurrence_matrix(df_aoi_occs_sub)

print("Computing co-occurence ...")
# Replace taxon keys by species name in the co-occurrence matrix
reclass_dict = dict(zip(species_keys, species_list))
df_aoi_cooc_f = df_aoi_cooc_f.rename(reclass_dict, axis='index').rename(reclass_dict, axis='columns')

# plot the percentage of co-occurrence of the first species in the species list to all other species in the species list
fig, axs = plt.subplots(1, 1, figsize=(8, 6))
n_counties_occ = df_aoi_cooc.loc[species_keys[0], species_keys[0]]
coocs = df_aoi_cooc_f[species_list[0]].sort_values(ascending=False)
# divide by total counties, in which target species is observed to get relative values
coocs = coocs/n_counties_occ * 100
coocs.drop(species_list[0]).plot(kind='bar')
plt.ylim(0, 100)
plt.ylabel("Percentage of co-occurrence with %s" % species_list[0])
plt.xticks(rotation=45)
plt.title(species_list[0])
plt.savefig(os.path.join(out_dir, species_list[0].replace(" ", "_") + "_co-occurrence.png"), bbox_inches='tight')
print("Saved figure to " + os.path.join(out_dir, species_list[0].replace(" ", "_") + "_co-occurrence.png"))
