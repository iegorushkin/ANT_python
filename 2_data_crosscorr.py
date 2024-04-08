# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 15:08:16 2021

@author: Igor
"""
import time
import os, sys
import numpy as np
import pandas as pd
import scipy.signal
from glob import glob
from obspy import read
from math import comb
import functions_crosscorr as f_cc
import functions_visualization as f_v

# Timer setup
# start = time.time()

# Define variables
# input_folder = 'E:/Work/Projects/Central_Kamchatka/MSEED/TEMP'
# input_folder = 'E:/Work/Projects/Central_Kamchatka/MSEED/PROCESSED'
# output_folder = 'E:/Work/Projects/Central_Kamchatka/Cross-correlations'
# pics_path = 'E:/Work/Projects/Central_Kamchatka/Cross-correlations/Pics/'

input_folder = 'H:/CK_PROCESSED_Z'
output_folder = 'H:/cross-correlations'
pics_path = 'H:/cross-correlations/pics/'

# Create a list of subfolders located in the input folder
subfolders_list = glob(input_folder + '/*')
# Number of subfolders
nsubfolders = len(subfolders_list)
# # Create a list of file lists in each subfolder
# subfolders_files = []
# Create a list with the number of files in each subfolder.
subfolders_nfiles = []
for subfolder in subfolders_list:
    temp = np.array(glob(subfolder + '/*'))
    # subfolders_files.append(temp)
    subfolders_nfiles.append(len(temp))
    del temp
# Create a list of ALL files located in the input folder.
# Format is: input folder -> station folder -> mseeds
files = np.array(glob(input_folder + '/*/*'))

# combinations_df
start1 = time.time()
# DataFrame with indices of pairs of stations
# If the list of file combinations for the required stations has already been
# evaluated, load it from the Excel file.
cc_combos_excel = 'H:/stations_combinations.xlsx'
if not os.path.exists(cc_combos_excel):
    cc_combos = f_cc.combinations_df_v3(files, nsubfolders,
                                        subfolders_nfiles, cc_combos_excel)
else:
    cc_combos = pd.read_excel(cc_combos_excel, index_col=0, header=(0, 1))
print('\nDuration of combinations_df_v3:')
end1 = time.time()
print(end1 - start1)

# Desired cross-correlation half length in seconds and samples
cc_len_sec = 300
cc_len = f_cc.corr_length_samples(cc_len_sec, files)

# Create DataFrame for containing cross-correlations averaged over all data
# 1) Determine how many stations are involved in the crosscorrelation
# nsubfolders
# 2) Count the total number of combinations of these stations
stations_ncombos = comb(nsubfolders, 2)
# Should be equal:
print('\nShould be equal:')
print(f"stations_ncombos - {stations_ncombos}"
      f"\ncc_combos.columns - {len(cc_combos.columns[::2])}")
# 3) Creating DataFrame
index = np.round(np.linspace(-cc_len_sec, cc_len_sec, 2*cc_len + 1), 1)
columns = [column[0] for column in cc_combos.columns[::2]]
cc_matrix = pd.DataFrame(np.zeros((2*cc_len + 1, stations_ncombos)),
                         index=index, columns=columns)

# Create array for containing hour-long cross-correlations
cc_matrix_hours = np.zeros((2*cc_len + 1, 24))

# Cross-correlation multiplication factor
factor = 100

# Flags
flag = 0
# %%
# Calculate and save daily cross-correlations

'''
cc_combos - DataFrame с мультииндексом столбцов.
Индекс верхнего уровня - название пары станций.
Ему соответствует два индекса нижнего уровня - название первой и второй станций
по-отдельности.
В столбцах, соответствующих индексам нижнего уровня содержатся, в свою очередь,
индексы строк из списка files.
Каждая строка files - путь файла mseed для загрузки и корреляции.
Таким образом, нужно:
1) рассмотреть каждый индекс верхнего уровня в cc_combos
2) поочередно извлечь каждую строку из 2 столбцов данных, соответствующих
индексам нижнего уровня
Каждая такая строка содержит индексы файлов из files, которые нужно
коррелировать.

Как делаем?

Двойной цикл:
- В первом делается проход по каждому второму значению cc_combos.columns.
Почему - можно понять из примера просто мультииндекса для двух станций:
MultiIndex([('CR09-CR10', 'CR09'),
            ('CR09-CR10', 'CR10')],
           )
Индекс верхнего уровня повторяется, нижнего нет.
Выполняя цикл таким образом, с помощью значения мультииндекса
cc_combos.columns[i-1] получаем доступ к столбцу для первой станции
i-той пары в cc_combos и при помощи cc_combos.columns[i] - для второй
станции i-той пары.
- во втором цикле поочередно извлекаются все пары значений cc_combos для
фиксированного индекса верхнего уровня, т.е. для фиксированной пары станций.
'''

for i in range(1, len(cc_combos.columns), 2):
    # Create DataFrame for containing day-long cross-correlations
    # Two-level MultiIndex for columns
    # same index as cc_matrix
    cc_matrix_days = pd.DataFrame(
        index=index, columns=pd.MultiIndex.from_product([[''], ['']]))

    for j in range(cc_combos.shape[0]):
        # Indices of files to load
        indx1 = int(cc_combos[cc_combos.columns[i-1]][j])
        indx2 = int(cc_combos[cc_combos.columns[i]][j])
        # Load data
        st1, sig1, st2, sig2 = f_cc.load_signals(files, indx1, indx2)

        # If both indices are equal to each other (as a result, equal to zero),
        # then there is no data for the current pair of stations and date.
        # Display a corresponding notification.
        if indx1 == indx2:
            print('\nNo data available!')
            print(f"Stations - {cc_combos.columns[i][0]}")
            print(f"Date - {cc_combos.index[j]}")
            print("Moving on to the next pair of indices...")
            # In order not to work with invalid / missing data:
            continue

        # Correlation
        else:
            print("\nNow cross-correlating... "
                  f"\nStations - {cc_combos.columns[i][0]},"
                  f" date - {cc_combos.index[j]}")
            # len(sig1) maybe != len(sig2), so mode='full' is used
            # (if mode = 'full', zero lag may not be included in the result)
            current_cc = scipy.signal.correlate(factor*sig1, factor*sig2,
                                                mode='full')
            # Index of zero-lag cross-correlation value
            # Zero lag - when zero value of sig1 is alligned
            # with zero value of sig2
            cc_zero = f_cc.corr_zerolag(sig1, sig2)

            # Check if the signal lengths are greater than
            # the cross-correlation window size
            if (len(sig1) < 2*cc_len + 1) or (len(sig2) < 2*cc_len + 1):
                print('The length of the signal (s) is smaller than'
                      ' the cross-correlation window size!')
                # Find a shorter signal and
                # set the correlation window length accordingly
                new_cc_len = min(len(sig1), len(sig2))
                # Check if cross-correlation is even or odd length,
                # trimm according to the its length
                if new_cc_len % 2 == 0:
                    new_cc_len = int(new_cc_len / 2)
                    current_cc = (
                        current_cc[(cc_zero - new_cc_len):(cc_zero + new_cc_len)]
                        )
                    # Add current_cc to the cc_matrix_hours
                    cc_matrix_hours[(cc_len - new_cc_len):(cc_len + new_cc_len),
                                    flag] = current_cc
                else:
                    new_cc_len = int(new_cc_len / 2)
                    current_cc = (
                        current_cc[(cc_zero - new_cc_len):(cc_zero + new_cc_len + 1)]
                        )
                    # Add current_cc to the cc_matrix_hours
                    cc_matrix_hours[(cc_len - new_cc_len):(cc_len + new_cc_len + 1),
                                    flag] = current_cc
                flag += 1
            # Case of longer signals
            else:
                # Trimm according to the desired length of cross-correlation
                current_cc = (
                    current_cc[(cc_zero - cc_len):(cc_zero + cc_len + 1)]
                    )
                # Add current_cc to the cc_matrix_hours
                cc_matrix_hours[:, flag] = current_cc
                flag += 1
        print('The current signals are cross-correlated')

        # Check if the next pair of files is from a different day.
        # Try to open the file following st1
        try:
            file_next = files[int(cc_combos[cc_combos.columns[i-1]][j+1])]
        except IndexError:
            print(f"\n***The end of the dataset is reached!***")
            # If we come to the end of the dataset,
            # then we average those hourly correlations that were collected,
            # even if there are not 24 of them.

            # Let's create a MultiIndex for the current temporary DataFrame
            # that will contain the average of the hourly cross-correlations
            # Station pair
            top_level = cc_combos.columns[i][0]
            # Date
            bottom_level = pd.to_datetime(str(st1[0].stats.starttime)[:10])
            bottom_level = bottom_level.strftime('%Y-%m-%d')  # trimm time
            # MultiIndex
            temp_columns = pd.MultiIndex.from_product([[top_level],
                                                       [bottom_level]])
            # Calculate the mean cross-correlation,
            # store it in the temporary DataFrame
            temp_df = pd.DataFrame(np.mean(cc_matrix_hours, axis=1),
                                   index=np.round(index, 1),
                                   columns=temp_columns)
            # Append it to cc_matrix_days_df
            cc_matrix_days = pd.concat([cc_matrix_days, temp_df], axis=1)

            # Reset some parameters
            cc_matrix_hours = np.zeros((2*cc_len + 1, 24))
            flag = 0
        else:
            # Not the end of the dataset
            # Collect month and day from current files.
            # Compare them with month and day from file_next.
            cur_md = str(st1[0].stats.starttime)[5:10]
            next_md = str(file_next[-14:-12] + '-' + file_next[-12:-10])
            if (cur_md != next_md):
                print("\n***The end of the current day is reached!***")
                # If true, then all files for the current day
                # are cross-correlated.

                # Let's create a MultiIndex for the current temporary DataFrame
                # that will contain the average of the hourly cross-correlation
                # Station pair
                top_level = cc_combos.columns[i][0]
                # Date
                bottom_level = (
                    pd.to_datetime(str(st1[0].stats.starttime)[:10]))
                bottom_level = bottom_level.strftime('%Y-%m-%d')  # trimm time
                # MultiIndex
                temp_columns = pd.MultiIndex.from_product([[top_level],
                                                           [bottom_level]])
                # Calculate the mean cross-correlation for the entire day,
                # store it in the temporary DataFrame
                temp_df = pd.DataFrame(np.mean(cc_matrix_hours, axis=1),
                                       index=np.round(index, 1),
                                       columns=temp_columns)
                # Append it to cc_matrix_days_df
                cc_matrix_days = pd.concat([cc_matrix_days, temp_df],
                                           axis=1)

                # Reset some parameters
                cc_matrix_hours = np.zeros((2*cc_len + 1, 24))
                flag = 0

    print("\nSaving daily cross-correlations"
          f" for the stations {cc_combos.columns[i][0]}...")
    # Save daily cross-correlations for the current station pair
    # Drop NaNs from the first column of this DF
    cc_matrix_days.dropna(axis=1, how='all', inplace=True)
    # ... And save in Excel format
    cc_matrix_days.to_excel(
        f"{output_folder}/{cc_combos.columns[i][0]}_daily.xlsx",
        float_format="%0.5f")
    print("Done!")
# %%
# Average daily cross-correlations and save the result
print('\nAveraging daily cross-correlations, saving the results...')

# Create a list of files with daily cross-correlations
excel_files = glob(output_folder + '/xlsx_daily/*-*')

# Loop through every file in f_ccile_list
for i in range(len(excel_files)):
    # Load current xlsx
    cc_matrix_days = pd.read_excel(excel_files[i], index_col=0,
                                   header=(0, 1))
    # Mean data from cc_matrix_days_df and save it to cc_matrix_df
    cc_matrix.iloc[:, i] = cc_matrix_days.groupby(level=0, axis=1).mean()

# Save cc_matrix_df - a matrix in the columns of which
# hourly cross-correlations averaged over the entire (available)
# period of operation of the corresponding pair of stations are stored
cc_matrix.to_excel(f"{output_folder}/averaged_crosscorrs_full.xlsx",
                   float_format="%0.5f")
print("Done!")
# %%
# Visualization
print("\nVisualization...")

# Check if cc_matrix_df is legit
# To do so, check if all values in this DataFrame are zeros
if (cc_matrix.iloc[:, :] == 0).all(axis=None):
    cc_matrix = pd.read_excel(
        f"{output_folder}/averaged_crosscorrs_full.xlsx", index_col=0)

# Create a list of files with daily cross-correlations
# f_ccile_list = glob(output_folder + '/*-*')

f_v.visualize_crosscorr(cc_matrix, excel_files, save_path=pics_path)
print("Done!")

# print('\nDuration of script execution:')
# end = time.time()
# print(end - start)
