#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
description
"""

__author__: "Hannah Weiser"
__email__: "h.weiser@stud.uni-heidelberg.de"

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from pygbif import occurrences as occ  # python client for the GBIF API
from pygbif import species
import sys
import os

plt.style.use('seaborn') # nicer figure style


def occs_to_gdf_clip(occ_dict, clip_gdf):
    occ_df = pd.DataFrame.from_dict(occ_dict['results'])
    occ_df = occ_df.dropna(subset=['decimalLatitude', 'decimalLongitude'])
    occ_df.to_csv("gbif_occ_" + occ_dict['results'][0]['scientificName'].replace(" ", "_") + ".csv")
    occ_df = occ_df.astype({'dateIdentified': 'datetime64', 'eventDate': 'datetime64'})
    occ_gdf = gpd.GeoDataFrame(occ_df, geometry=gpd.points_from_xy(occ_df.decimalLongitude, occ_df.decimalLatitude))
    clip_crs = clip_gdf.crs
    occ_gdf.set_crs(clip_crs, inplace=True)
    occ_gdf_clipped = gpd.clip(occ_gdf, clip_gdf)

    return occ_gdf_clipped


def occs_to_gdf(species_dict):
    occ_df = pd.DataFrame.from_dict(species_dict['results'])
    occ_df = occ_df.dropna(subset=['decimalLatitude', 'decimalLongitude'])
    #occ_df = occ_df.astype({'dateIdentified': 'datetime64', 'eventDate': 'datetime64'})
    occ_gdf = gpd.GeoDataFrame(occ_df, geometry=gpd.points_from_xy(occ_df.decimalLongitude, occ_df.decimalLatitude))

    return occ_gdf


def sjoin_and_merge(polys, points, pivot_index, count_column):
    crs = polys.crs
    points.set_crs(crs, inplace=True)
    dfsjoin = gpd.sjoin(polys, points)
    dfpivot = pd.pivot_table(dfsjoin, index=pivot_index, columns=count_column, aggfunc={count_column: len})
    dfpivot.columns = dfpivot.columns.droplevel()
    polys_new = polys.merge(dfpivot, how='left', on=pivot_index)

    return polys_new


# input data
species_list = ['Sphex funerarius',
                'Bicolorana bicolor',
                'Conocephalus fuscus',
                'Conocephalus dorsalis',
                'Decticus verrucivorus',
                'Gryllus campestris',
                'Meconema meridionale',
                'Meconema thalassinum',
                'Metrioptera roeselii',
                'Phaneroptera nana',
                'Phaneroptera falcata',
                'Platycleis albopuntata',
                'Tesselana tesselata',
                'Tettigonia viridissima'
                ] #read from file

fp = r"./NUTS_RG_10M_2021_4326_LEVL_3.shp/NUTS_RG_10M_2021_4326_LEVL_3.shp" #read from file
out_dir = "./output"

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

gdf_aoi_adm = gpd.read_file(fp)
gdf_aoi_adm = gdf_aoi_adm[gdf_aoi_adm['CNTR_CODE']=='DE']
gdf_aoi_occs = gdf_aoi_adm.copy()

for sp in species_list:
    occ_dict = occ.search(scientificName=sp, country="DE")
    if occ_dict['results'] == []:
        continue
    occ_gdf = occs_to_gdf(occ_dict)
    occ_gdf.to_csv("gbif_occ_" + occ_dict['results'][0]['scientificName'].replace(" ", "_") + ".csv")
    occ_gdf.plot(column="year")
    gdf_aoi_occs = sjoin_and_merge(gdf_aoi_occs, occ_gdf, 'NUTS_ID', 'taxonKey')

species_keys = [species.name_backbone(x)['usageKey'] for x in species_list]

n_plots = len([key for key in species_keys if key in gdf_aoi_occs.columns])
n_cols = 4
n_rows = int(np.ceil(n_plots/n_cols))
fig, axs = plt.subplots(n_rows,n_cols, figsize=(20,10))
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

