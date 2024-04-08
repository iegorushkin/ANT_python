# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 17:53:05 2021

@author: Igor
"""
import time
import os
import sys
import functions_visualization as f_v
import functions_filtering as f_f
from glob import glob
from obspy import read, read_inventory
sys.path.append(r'E:\Work\scripts_python\python_external_packages\greentools-master')
from greentools.core import write_st_to_mseed

'''
This script:
-  pre-filters, removes the instrument response and mean from the records
of stations;
- demeans and detrends again;
- filters the signal with a bandpass Butterworth filter;
- decreases the sampling rate;
- normalizes data in frequency and time domains
- saves the result in mseed format
'''

# Define variables

# input_folder = r'E:/Work/Projects/Central_Kamchatka/MSEED/RAW/1hour'
# output_folder = r'E:/Work/Projects/Central_Kamchatka/MSEED/PROCESSED'
# pics_path = r'E:/Work/Projects/Central_Kamchatka/data_processing_pics/'
input_folder = 'H:/CK_RAW_Z/1hour/perm_stations'
# input_folder = 'H:/temp'
output_folder = 'H:/CK_PROCESSED_Z'
pics_path = 'H:/data_processing_pics/'
inventory_path = 'H:/Central_Kamchatka/Meta_data/inventory_perm_stations.xml'
start = time.time()

goal_sampling_rate = 10

# Permanent stations
pre_filt_bh = (0.018, 0.022, 1, 2)
# bh_list = ['APC', 'AVH', 'DAL', 'IVS', 'KRM', 'PET', 'SBLV', 'SPN', 'UBL']
bh_list = ['APC', 'DAL', 'IVS', 'KRM', 'PET', 'SBLV', 'SPN', 'UBL']
# Temporary stations
pre_filt_mc = (0.005, 0.01, 1, 2)
mc_list = ['CR05', 'CR06', 'CR07', 'CR08', 'CR09', 'CR10', 'CR11', 'CR12',
           'CR13', 'CR14', 'WH03', 'WR02']
pre_filt_guralp = (0.030, 0.035, 1, 2)
guralp_list = ['CR01', 'CR02', 'CR03', 'CR04', 'EH01', 'EH02', 'EH04', 'EH05',
               'EH06', 'EH07', 'EH08', 'EH09', 'EH10', 'EH11', 'WH02', 'WH04',
               'WR01', 'WR03', 'WR04', 'WR05']

# Frequency band of interest
frange = [0.035, 1.0]

# # Clean folder with pics if necessary
# for f in os.listdir(pics_path):
#     os.remove(os.path.join(pics_path, f))

# Reading inventory, creating a list of files located in the input folder.
inv = read_inventory(inventory_path)
file_list = glob(input_folder + '/*/*')  # \ instead of /
# print(file_list)

# Flag for visualization of the processing steps
flag = 1000
flag_divider = 10000
# %%
# Calculations

# Cycling through every file in file_list
for file in file_list:
    # Saving names of the current station and the current file
    station, file_name = file.split('\\')[-2:]

    # Reading the current file
    st = read(file)

    # Checking what type of sensor was used in this station
    if st[0].stats.station in bh_list:
        response_prefilter = pre_filt_bh
    elif (st[0].stats.station in mc_list) or (st[0].stats.station in guralp_list):
        print("The current data is taken from a temporary station. Skip...")
        continue
    else:
        print(f"Not sure what type of sensor is used in {st[0].stats.station}")
        continue
        # sys.exit()

    # Visualize raw st
    if flag % flag_divider == 0:
        # Plot of the original signal
        f_v.visualize_signal(st, freq_xlim=[0, 2], dont_open=True,
                             file_name=file_name, save_path=pics_path)

    # Pre-filtering, removing instrument response, mean
    # (before and after deconvolution) and linear trend
    remove_response_pic = (pics_path + file_name + '_remove_response.png')
    if flag % flag_divider == 0:
        st.remove_response(inventory=inv, pre_filt=response_prefilter,
                           water_level=None, zero_mean=True,
                           taper=True, taper_fraction=0.04,
                           plot=remove_response_pic)
    else:
        st.remove_response(inventory=inv, pre_filt=response_prefilter,
                           water_level=None, zero_mean=True,
                           taper=True, taper_fraction=0.04,)
    if flag % flag_divider == 0:
        # Plots the signal after remove_response;
        # also shows this signal without mean and linear trend.
        f_v.visualize_signal_adv(st, freq_xlim=[0, 2], dont_open=True,
                                 file_name=file_name, save_path=pics_path)
    # st.detrend uses scipy.signal.detrend, so this is basically the same thing
    # that is calculated in vis_sig
    st.detrend('demean')
    st.detrend('linear')

    # Applying lowpass and highpass filters
    st = f_f.passband_butter(st, frange, 4, 9)
    # Plot the result of filtration
    if flag % flag_divider == 0:
        f_v.visualize_fsignal(st, frange, 4, 9, freq_xlim=[0, 2],
                              dont_open=True, file_name=file_name,
                              save_path=pics_path)

    # Resampling to the goal sampling rate
    decimate_factor = int(st[0].stats.sampling_rate / goal_sampling_rate)
    st.decimate(decimate_factor, no_filter=True)
    # Plot the result of downsampling
    if flag % flag_divider == 0:
        f_v.visualize_signal(st, freq_xlim=[0, 2], dont_open=True,
                             file_name=file_name, save_path=pics_path)

    ## Uncomment for 1 bit normalization
    # Spectral normalization
    # st = f_f.spect_whiten(st, frange, file_name)
    # # Plot
    # if flag % flag_divider == 0:
    #     f_v.visualize_signal(st, freq_xlim=[0, 2], dont_open=True,
    #                          file_name=file_name, save_path=pics_path)
    ##

    # Time domain normalization
    # st = f_f.onebit_norm(st)
    st = f_f.ram_norm(st, 1/frange[0])
    # Plot
    if flag % flag_divider == 0:
        f_v.visualize_signal(st, freq_xlim=[0, 2], dont_open=True,
                             file_name=file_name, save_path=pics_path,)
                             # time_xlim=[0, 1])

    # Spectral normalization
    st = f_f.spect_whiten(st, frange, file_name)
    # Plot
    if flag % flag_divider == 0:
        f_v.visualize_signal(st, freq_xlim=[0, 2], dont_open=True,
                             file_name=file_name, save_path=pics_path)
    flag += 1

    # Saving result in a mseed file as 32-bit floats - this is the same
    # as SAC precision
    output_file = os.path.join(output_folder, station, file_name)
    write_st_to_mseed(st, output_file)

print('\nDuration of script execution:')
end = time.time()
print((end - start)/60)
