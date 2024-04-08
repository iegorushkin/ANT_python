# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 17:16:27 2021

@author: Igor
"""
import numpy as np
import pandas as pd
from obspy.geodetics import gps2dist_azimuth

# Load a file with cross-correlations
crosscorrs_path = 'H:/cross-correlations/averaged_crosscorrs_full.xlsx'
crosscorrs = pd.read_excel(crosscorrs_path, index_col=0)
# Create an array with the names of all pairs of stations and count them
station_pairs = crosscorrs.columns.values
nstation_pairs = len(station_pairs)
#  Create an array to store the calculated distances between stations
dist_array = np.zeros(nstation_pairs)

# Load the file with the location of each station
coordinates_path = 'H:/stations_locations.xlsx'
coordinates = pd.read_excel(coordinates_path, index_col=0)

for i in range(nstation_pairs):
    cur_stat1, cur_stat2 = station_pairs[i].split('-')

    # Find latitude and longitude of cur_stat1 in coordinates
    lat1 = coordinates['latitude'][cur_stat1]
    lon1 = coordinates['longitude'][cur_stat1]
    # Find latitude and longitude of cur_stat2 in coordinates
    lat2 = coordinates['latitude'][cur_stat2]
    lon2 = coordinates['longitude'][cur_stat2]

    # Calculate the distance between cur_stat1 and cur_stat2
    dist_array[i], az, baz = gps2dist_azimuth(lat1, lon1, lat2, lon2)

# Combine the results in a DataFrame and save it as an .xlsx document
pairs_dist = pd.DataFrame(data=np.round(dist_array).astype(int),
                          index=station_pairs, columns=['Distance, m'])
pairs_dist.to_excel('H:/stations_distance.xlsx')
