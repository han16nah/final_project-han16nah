#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
description
"""

__author__: "Hannah Weiser"
__email__: "h.weiser@stud.uni-heidelberg.de"
__date__: "11/2020"

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


def get_nuts(spatialtype='RG', resolution='10M', year='2021', projection='4326', subset='LEVL_3', output_file = "NUTS"):
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
    format = 'geojson'
    file = "%s_%s_%s_%s_%s_%s.%s" % (theme, spatialtype, resolution, year, projection, subset, format)
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


def occs_to_gdf(occurence_dict):
    """
    Converts an occurence dictionary to a geopandas GeoDataFrame and projects it to WGS84 (crs of occurences in GBIF).
    :param occurence_dict: occurence dictionary as returned by the GBIF API (occurence section)
    :return: Occurence GeoDataFrame
    """
    occ_df = pd.DataFrame.from_dict(occurence_dict['results'])
    occ_df = occ_df.dropna(subset=['decimalLatitude', 'decimalLongitude'])
    occ_df = occ_df.astype({'year':'int'})
    occ_gdf = gpd.GeoDataFrame(occ_df, geometry=gpd.points_from_xy(occ_df.decimalLongitude, occ_df.decimalLatitude))
    occ_gdf.set_crs('EPSG:4326', inplace=True)

    return occ_gdf


def clip_gdf(gdf, clip_feat, crs = None):
    """
    Clips a dataframe by a clip feature and makes sure they are both in the same crs.
    :param gdf: Input GeoDataFrame to clip.
    :param clip_feat: Clip feature (GeoDataFrame or GeoSeries)
    :param crs: Coordinate reference system to project the input data to, in the form "EPSG:4326"
        If "None", will use the crs of gdf
    :return: Clipped GeoDataFrame
    """
    # Catch the error if one of the input DataFrames has no coordinate reference system
    try:
        if crs == None:
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


def build_cooccurence_matrix(df):
    """
    Reclassifies the values in a dataframe and then builds  a co-occurence matrix.
    Values > 0      are reclassified to     1
    NaN-values      are reclassified to     0
    :param df: DataFrame to build the co-occurence matrix for
    :return: Tuple:
        (Co-occurence matrix, Co-occurence matrix with diagonals filled with 0)
    """
    # reclassify: values > 0 --> 1; NaN --> 0
    df.values[df.values > 0] = 1
    df.fillna(0, inplace=True)
    df_asint = df.astype(int)
    df_cooc = df_asint.T.dot(df_asint)
    df_cooc_f = df_cooc.copy()
    np.fill_diagonal(df_cooc_f.values, 0)

    return df_cooc, df_cooc_f

# input data
config_file = sys.argv[1]
with open(config_file) as src:
    config = json.load(src)

species_list = config["species_list"]
out_dir = config["output_dir"]
try:
    country_code = config["country_code"]
except:
    country_code = None

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

gdf_aoi_adm = get_nuts()
if country_code != None:
    gdf_aoi_adm = gdf_aoi_adm[gdf_aoi_adm['CNTR_CODE'] == country_code]
gdf_aoi_occs = gdf_aoi_adm.copy()

species_keys = [species.name_backbone(x)['usageKey'] for x in species_list]

i = 0
for sp in species_list:
    occ_dict = occ.search(scientificName=sp, country=country_code)
    if occ_dict['results'] == []:
        continue
    occ_gdf = occs_to_gdf(occ_dict)

    occ_gdf.to_csv("gbif_occ_" + occ_dict['results'][0]['scientificName'].replace(" ", "_") + ".csv")
    fig, ax = plt.subplots(figsize=(10,10))
    gdf_aoi_adm.plot(ax=ax, color='lightgrey')
    gdf_aoi_adm.geometry.boundary.plot(ax=ax, alpha=0.25, color='black', lw=0.5)
    occ_gdf.plot(column="year", ax=ax, cmap="OrRd", legend=True, categorical=True)
    plt.title("%s - Observations by year" % sp)
    plt.savefig(os.path.join(out_dir, sp + "_occ_by_year.png"))
    plt.clf()
    gdf_aoi_occs = sjoin_and_merge(gdf_aoi_occs, occ_gdf, 'NUTS_ID', 'taxonKey')

n_plots = len([key for key in species_keys if key in gdf_aoi_occs.columns])
n_cols = 4
n_rows = int(np.ceil(n_plots/n_cols))
fig, axs = plt.subplots(n_rows,n_cols, figsize=(20,15))
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
plt.clf()

keys = ['NUTS_ID'] + species_keys
keys = [key for key in keys if key in gdf_aoi_occs.columns]
df_aoi_occs_sub = gdf_aoi_occs[keys]
df_aoi_occs_sub.set_index('NUTS_ID', inplace=True)
df_aoi_cooc, df_aoi_cooc_f = build_cooccurence_matrix(df_aoi_occs_sub)

reclass_dict = dict(zip(species_keys, species_list))
df_aoi_cooc_f = df_aoi_cooc_f.rename(reclass_dict, axis='index').rename(reclass_dict, axis='columns')

fig, axs = plt.subplots(1,1, figsize=(8,6))
n_counties_occ = df_aoi_cooc.loc[species_keys[0], species_keys[0]]
coocs = df_aoi_cooc_f[species_list[0]].sort_values(ascending=False)
coocs = coocs/n_counties_occ * 100
coocs.drop(species_list[0]).plot(kind='bar')
plt.ylim(0, 100)
plt.ylabel("Percentage of co-occurence with %s" % species_list[0])
plt.xticks(rotation=45)
plt.title(species_list[0])
plt.savefig(os.path.join(out_dir, species_list[0].replace(" ", "_") + "_co-occurence.png"), bbox_inches='tight')
