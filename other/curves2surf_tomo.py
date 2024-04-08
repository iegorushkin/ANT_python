# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 13:18:04 2022

@author: Igor
"""
# Imports and setting-up basic variables
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.interpolate
from glob import glob



input_folder = 'D:\\Central_Kamchatka\\dispersion curves\\raw'
file_list = glob(input_folder + '\\*')
file_list.sort()

output_folder = 'D:\\Central_Kamchatka\\dispersion curves\\'
# %%
# Load and assemble all dispersion curves in one DataFrame object
df = pd.DataFrame(columns=pd.MultiIndex.from_product([[''], ['']]))
for file in file_list:
    cur_df = pd.read_csv(file, sep="\t", header=None)
    cur_stations = file.split(sep='\\')[-1].split(sep='_')[0]
    cur_df.columns = pd.MultiIndex.from_product([[cur_stations],
                                                 ['period', 'group velocity']])
    df = pd.concat([df, cur_df], axis=1)
df.dropna(axis=1, how='all', inplace=True)

# Load and save station information
stations_info = pd.read_excel('D:\\Central_Kamchatka\\stations_locations.xlsx',
                              index_col=0)
stations_info['elevation'] = np.round(stations_info['elevation'], 4)
# %%
# Interpolation of dispersion curves
idf = pd.DataFrame(data=np.zeros((100, df.shape[1])), columns=df.columns)

# Loop through every top level of df MultiIndex
for i in range(0, df.shape[1], 2):
# for i in range(0, 4, 2):
    x = df.iloc[:, i].dropna()
    y = df.iloc[:, i+1].dropna()

    # Get extrapolation classes
    f_lin = scipy.interpolate.interp1d(x.values, y.values, kind='linear',
                                       fill_value='extrapolate')

    # Set the first and last point of extrapolated data.
    # The coordinates of the first point are rounded down to the nearest int,
    # and the coordinates of the last point are rounded up to the nearest int.
    # fpoint = np.floor(x.values[0] / 0.5) * 0.5
    # lpoint = np.ceil(x.values[-1] / 0.5) * 0.5
    fpoint = np.floor(x.values[0])
    lpoint = np.ceil(x.values[-1])

    # Extrapolate!
    # x_lin = np.arange(fpoint, lpoint + 0.5, 0.5)
    x_lin = np.arange(fpoint, lpoint+1)
    y_lin = f_lin(x_lin)
    length = len(x_lin)

    # Save to DataFrame
    # Some of the first values look like artifacts,
    # this needs to be taken into account.
    if y_lin[0] > 4 or x_lin[0] == 0:
        idf.iloc[:length-1, i] = x_lin[1:]
        idf.iloc[:length-1, i+1] = y_lin[1:]
    else:
        idf.iloc[:length, i] = x_lin
        idf.iloc[:length, i+1] = y_lin

idf = idf[(idf != 0).any(axis=1)]  # remove unneeded zeros

# temp = idf.xs(key='group velocity', axis=1, level=1)
# temp_max = temp.max(axis=0)
# bad_result = pd.concat([df['WR01-WR02'], idf['WR01-WR02']], axis=1)
# %%
# Plotting
fig, ax = plt.subplots(2, 1, sharex=True,
                       gridspec_kw={'height_ratios': [2, 1]}, figsize=(12, 6))
# %%
# Let's plot all those curves!

# Create DF where indices are periods, columns are group velocities
df_final = pd.DataFrame(data=[], index=range(1, 41))
df_final.index.name = 'Periods'
df_final.reset_index(inplace=True)

# Loop through every top level of df MultiIndex
for i in range(0, df.shape[1], 2):
    # First, need to remove zeros from the end of the temp DataFrame
    temp_df = idf.iloc[:, i:i+2]
    temp_df = temp_df[(temp_df != 0).any(axis=1)]
    # And plot
    ax[0].plot(temp_df.iloc[:, 0], temp_df.iloc[:, 1], 'k', alpha=0.33)
    # Make plot pretty
    ax[0].set(xlim=[1, 35], ylim=[0, 5])
    # ax[0].set_ylabel('Group velocity, km/s', fontsize=12)
    # ax[0].set_title('All manually picked dispersion curves', fontsize=14)
    # for RNF visualization
    ax[0].set_ylabel('групповая скорость, км/с', fontsize=12)

    ## This block is for use in further cells ##
    # Preparing temp_df for merging...
    temp = temp_df.droplevel(0, 1).rename(
        columns={'group velocity': temp_df.columns[0][0],
                 'period': 'Periods'})
    temp.iloc[:, 0].astype(int)
    # Merging!
    df_final = pd.merge(df_final, temp, on='Periods', how='left')
    df_final.set_index('Periods', inplace=True)

ax[0].grid()
fig.tight_layout()
# %%
# Now, counts and plot how many curves there are at a specific period
# period_count = pd.Series(data=np.zeros(40), index=range(1, 41))
period_count = pd.DataFrame(data=np.zeros(40), columns=['number of occurrences'],
                            index=range(1, 41))
for i in range(1, 41):
    period_count.loc[i] = np.sum(idf.values == i)

# Save period_count
period_count.to_csv('D:\\Central_Kamchatka\\period_count.txt', sep='\t')

ax[1].plot(period_count.index, period_count.values, '*k', markersize=10)
# ax[1].set_title('Number of observations per period', fontsize=14)
# ax[1].set_ylabel('N', fontsize=12)
# ax[1].set_xlabel('Periods, s', fontsize=12)
ax[1].set(ylim=[0, 515])
ax[1].set_ylabel('N', fontsize=12)
ax[1].set_xlabel('периоды, с', fontsize=12)
ax[1].grid()
fig.tight_layout()
#%%
# Save the figure
save_path = 'E:/Work/Projects/Central_Kamchatka/_Paper/Figures/'
output_file = save_path + "rays_stats.png"
fig.savefig(output_file, bbox_inches='tight', dpi=300)
# %%
# Let's determine with what periods the tomographic 'actions' will be performed
period_count['% compared to the most frequent period'] = (
    period_count.iloc[:, 0] / period_count.iloc[:, 0].max())
# Throw out periods that are 6.66 times less common than the most common one
periods_final = (
    period_count[period_count['% compared to the most frequent period'] >= 0.15])
# %%
del x, x_lin, y, y_lin, cur_df, cur_stations, fpoint, lpoint, length, temp_df

# Loop throught every relevant period
for i in range(1, 31):
    # Extract existing group velocities for period i
    cur_ser = df_final.loc[i].dropna()
    # print(temp_ser)

    # Create an empty matrix to store the coordinates of the first and second
    # stations and the group velocity between them on period i.
    cur_rays = np.zeros((cur_ser.shape[0], 7))

    # Loop through all pairs of stations stored in cur_ser
    for idx, stations in enumerate(cur_ser.index.values):
        # 'Split' the current pair of stations
        station_a = stations.split(sep='-')[0]
        station_b = stations.split(sep='-')[1]

        # coordinates of station 1
        cur_rays[idx, :3] = stations_info.loc[station_a].values
        # coordinates of station 1
        cur_rays[idx, 3:6] = stations_info.loc[station_b].values
        # group velocity between them
        cur_rays[idx, 6] = cur_ser[stations]

    # Saving cur_rays
    if i < 10:
        np.savetxt(fname=output_folder + '\\inidata\\rays0' + str(i) + '.dat',
                   X=cur_rays, fmt='%.4f', delimiter='\t')
    else:
        np.savetxt(fname=output_folder + '\\inidata\\rays' + str(i) + '.dat',
                   X=cur_rays, fmt='%.4f', delimiter='\t')

np.savetxt(fname=output_folder + '\\inidata\\periods.dat',
           X=np.arange(1, 31), fmt='%.0f')
