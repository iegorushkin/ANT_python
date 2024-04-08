# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 16:48:41 2021

@author: Igor
"""
import pandas as pd
import numpy as np
import os
import functions_crosscorr as f_cc
# %%
# Different inputs
input_folder = 'H:/CK_PROCESSED_Z'
output_folder = 'H:/cross-correlations/SAC'

# load DataFrames
cc_matrix = pd.read_excel("H:/cross-correlations/averaged_crosscorrs_full.xlsx",
                          index_col=0)
stations_info = pd.read_excel(input_folder + "/stations_locations.xlsx",
                              index_col=0)
# Check DataFrames
print(f"\n{cc_matrix.head()}")
print(f"\n{stations_info.head()}")

# Create appropriate DataFrames for the negative and positive portions of
# the cross-correlations; also for the average
new_len = int(np.ceil(cc_matrix.shape[0]/2))
cc_matrix_n = cc_matrix.iloc[new_len-1::-1]  # .reset_index(drop=True)
cc_matrix_p = cc_matrix.iloc[new_len-1:]
cc_matrix_avg = pd.concat([cc_matrix_n.reset_index(drop=True),
                           cc_matrix_p.reset_index(drop=True)], axis=1)
cc_matrix_avg = cc_matrix_avg.groupby(axis=1, level=0).mean()
# same as
# cc_matrix_avg1=cc_matrix_avg.groupby(by=cc_matrix_avg.columns, axis=1).mean()
# (cc_matrix_avg == cc_matrix_avg2).all(axis=0)
cc_matrix_avg.set_index(cc_matrix_p.index, inplace=True)  # баловство, но пусть

# Create a dictionary whose keys and values correspond
# to different types of output
types = {'negative': cc_matrix_n, 'positive': cc_matrix_p,
         'average': cc_matrix_avg}
# %%
# Create a dictionary where keys are station names
# and values are their obspy Stats
stats_dict = f_cc.create_stats_dict(input_folder)

# Loop through keys and values in types dictionary
for cur_type, cur_df in types.items():
    # Loop through every column in the cross-correlation matrix
    # In each loop split the column name and use its parts
    # to create an ObSpy Trace object
    for column in cur_df.columns:
        # print("\n" + column)
        # print(column.split('-'))
        station1, station2 = column.split('-')

        # Creating Obspy Trace object based on various sources
        cur_tr = f_cc.corr_data_to_trace(cur_df[column].to_numpy(),
                                         station1, station2, stations_info,
                                         stats_dict, sac=True)

        # Save this object
        # Diffirent type of the output structure
        # filename = f"{output_folder}/SAC/{cur_type}/{column}_{cur_type}.sac"
        # cur_tr.write(filename, format="SAC")

        # Check if a directory for the current output is created
        # Create it if it does not
        cur_output_dir = os.path.join(output_folder, column)
        if (not os.path.isdir(cur_output_dir)):
            os.mkdir(cur_output_dir)
        # Saving the result...
        filename = cur_output_dir + f"/{column}_{cur_type}.sac"
        cur_tr.write(filename, format="SAC")
#%%
# from obspy import read

### TESTING ###
# test_sac1 = read("Cross-correlations/average/CR13-EH09_average.sac")
# test_sac2 = read("Cross-correlations/negative/CR13-EH09_negative.sac")
# test_sac3 = read("Cross-correlations/positive/CR13-EH09_positive.sac")

# print('')
# print(test_sac1[0].stats)
# print(test_sac2[0].stats)
# print(test_sac3[0].stats)

# test_mseed = read("E:/Work/Projects/Central_Kamchatka/test_mseed")
# print(test_mseed[0].stats)

# filename2 = column + '_type.mseed'
# cur_tr.write(filename2, format="mseed")
