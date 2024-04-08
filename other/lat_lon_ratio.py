# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 12:23:23 2022

@author: Igor
"""
import numpy as np

# mean latitude of the study area (in degrees)
mlat = 53.5
# mean latitude of the study area (in radians)
mlat_d = mlat*np.pi/180

# Source: https://en.wikipedia.org/wiki/Geographic_coordinate_system
# WGS84 spheroid
# calculate the length in meters of a degree of latitude at latitude mlat
len_lat = (111132.92 - 559.82*np.cos(2*mlat_d)
           + 1.175*np.cos(4*mlat_d) - 0.0023*np.cos(6*mlat_d))
# calculate the length in meters of a degree of longitude at latitude mlat
len_lon = (111412.84*np.cos(mlat_d) - 93.5*np.cos(3*mlat_d)
           + 0.118*np.cos(5*mlat_d))
# calculate ratio between len_lat and len_lon
ratio = np.round(len_lat / len_lon, 3)

print(f"The length in meters of a degree of latitude at latitude {mlat}")
print(np.round(len_lat/1000, 3))

print(f"\nThe length in meters of a degree of longitude at latitude {mlat}")
print(np.round(len_lon/1000, 3))

print('\nRatio: {}'.format(ratio))






# # mean latitude and longitude of the study area (in degrees)
# mlat = 53.5
# mlon = 53.5
# # mean latitude and longitude of the study area (in radians)
# mlat_d = mlat*np.pi/180
# mlon_d = mlon*np.pi/180

# # Source: https://en.wikipedia.org/wiki/Geographic_coordinate_system

# # calculate the length in meters of a degree of latitude at latitude mlat
# len_lat = (111132.92 - 559.82*np.cos(2*mlat_d)
#            + 1.175*np.cos(4*mlat_d) - 0.0023*np.cos(6*mlat_d))
# # calculate the length in meters of a degree of longitude at longitude mlat
# len_lon = (111412.84*np.cos(mlon_d) - 93.5*np.cos(3*mlon_d)
#            + 0.118*np.cos(5*mlon_d))
# # calculate ratio between len_lat and len_lon
# ratio = np.round(len_lat / len_lon, 3)

# print(f"The length in meters of a degree of latitude at latitude {mlat}")
# print(np.round(len_lat/1000, 3))

# print(f"\nThe length in meters of a degree of longitude at longitude {mlon}")
# print(np.round(len_lon/1000, 3))

# print('\nRatio: {}'.format(ratio))
