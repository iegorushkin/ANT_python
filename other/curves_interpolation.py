# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 22:15:21 2021

@author: Igor
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.interpolate
from glob import glob
# %% Testing different kinds of interpolation

# Load data
data = np.loadtxt("H:/dispersion curves/raw/PET-UBL_average__VGDC.txt")

# Plot it
fig, ax = plt.subplots(figsize=(8, 6))
ax.set(xlabel='Period, sec', ylabel='Group velocity, km/s',
       title='Example of a dispersion curve')
ax.grid()

# Get extrapolation classes
f_lin = scipy.interpolate.interp1d(data[:, 0], data[:, 1], kind='linear',
                                   fill_value='extrapolate')
f_qd = scipy.interpolate.interp1d(data[:, 0], data[:, 1], kind='quadratic',
                                  fill_value='extrapolate')
f_cub = scipy.interpolate.interp1d(data[:, 0], data[:, 1], kind='cubic',
                                   fill_value='extrapolate')

# Set the first and last point of extrapolated data.
# The coordinates of the first point are rounded down to the nearest .5,
# and the coordinates of the last point are rounded up to the nearest .5.
fpoint = np.floor(data[0, 0] / 0.5) * 0.5
lpoint = np.ceil(data[-1, 0] / 0.5) * 0.5

# Extrapolate!
x = np.arange(fpoint, lpoint + 0.5, 0.5)
y_lin = f_lin(x)
y_qd = f_qd(x)
y_cub = f_cub(x)

ax.plot(x, y_lin, '*-', label='linear', lw=3, markersize=7)
ax.plot(x, y_qd, 'x-', label='quadratic', lw=3, markersize=7)
ax.plot(x, y_cub, '^-', label='cubic', lw=3, markersize=7)
ax.plot(data[:, 0], data[:, 1], 's-', label='raw', lw=3, markersize=7)
ax.legend()
# %% Let's use cubic!
folder_list = glob('Cross-correlations/Curves/' + '/*/')
exclude = folder_list[2]  # pics
fnumber = 15  # Number of curves in every folder

file_lists = []  # list of file lists from every folder (avg, neg, pos)
for i in range(len(folder_list)):
    if folder_list[i] != exclude:
        cur_file_list = glob(folder_list[i] + '*')
        file_lists.append(cur_file_list)

# Loop through the number of curves in each folder
for i in range(fnumber):
    print()
    # Prepare figure
    cur_station_names = file_lists[0][i][-27:-18]
    fig, ax = plt.subplots(figsize=(12, 4))
    fig.suptitle("Group velocity dispersion curves inferred from stations "
                 + f"{cur_station_names}", fontsize=15, y=0.94)
    plt.close(fig)  # do not display figure
    ax.set_xlabel('Period, s', fontsize=12)
    ax.set_ylabel('Group velocity, km/s', fontsize=12)
    ax.grid(which='both', linewidth=0.5, color='k')
    # Create list of labels and colors for future use
    llabels = ['average', 'negative', 'positive', ]
    lcolors = ['magenta', 'red', 'blue', ]

    # For each pair of stations, several types of curves are set,
    # which are saved in separate folders. Loop through each folder.
    for j in range(len(file_lists)):
        print(file_lists[j][i])
        cur_data = np.loadtxt(file_lists[j][i])
        ax.plot(cur_data[:, 0], cur_data[:, 1], c=lcolors[j],
                label=llabels[j], lw=2)
    ax.legend(title="       Used part of\nthe cross-correlation:", fontsize=None,
              framealpha=1)
    # Save figure
    fig_filename = ('Cross-correlations/Curves/pics/'
                    + cur_station_names + '.png')
    fig.savefig(fig_filename, bbox_inches='tight')
