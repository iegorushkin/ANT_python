# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 15:50:36 2021

@author: Igor
"""
import os
from obspy import read, UTCDateTime

# input_dir = r'H:\1day'
# output_dir = r'H:\1hour'

input_dir = r'H:\CK_RAW_Z\1day\perm_stations'
output_dir = r'H:\CK_RAW_Z\1hour\perm_stations'

# List of folders in input directory
folder_list = os.listdir(input_dir)
# folder_list = ['AVH_BH'] - no response!

# Go through every folder in input directory
for folder in folder_list:
    folder_path = os.path.join(input_dir, folder)
    # List of files in current folder
    infilist = os.listdir(folder_path)
    # Go through every file in current folder
    for input_file in infilist:
        # Create template for reading the current MSEED-file
        st = read()
        st.clear()
        # Define full path to the current input file
        input_full_name = os.path.join(folder_path, input_file)
        print(f'\nFull path to the current input file: \n{input_full_name}')
        # Reading the current file
        st = read(input_full_name)

        # Now st will be cut into 1 hour-long mseed files
        # Define start time and end time of st
        starttime = UTCDateTime(year=st[0].stats.starttime.year,
                                julday=st[0].stats.starttime.julday,
                                hour=st[0].stats.starttime.hour,
                                minute=st[0].stats.starttime.minute,
                                second=st[0].stats.starttime.second
                                )
        print(f'\nStart time:\n{starttime}')
        endtime = UTCDateTime(year=st[0].stats.endtime.year,
                              julday=st[0].stats.endtime.julday,
                              hour=st[0].stats.endtime.hour,
                              minute=st[0].stats.endtime.minute,
                              second=st[0].stats.endtime.second
                              )
        print(f'\nEnd time:\n{endtime}')

        # Let's cut!
        flag = True
        cur_starttime = starttime
        while flag:
            cur_endtime = UTCDateTime(year=cur_starttime.year,
                                      julday=cur_starttime.julday,
                                      hour=cur_starttime.hour,
                                      minute=59, second=59, microsecond=99*10**4)
            print(f'\nCurrent start time is {cur_starttime}')
            print(f'Current end time is {cur_endtime}')

            # Cycle through every trace that will be saved
            for tr in st:
                tr_sliced = tr.slice(starttime=cur_starttime,
                                     endtime=cur_endtime,
                                     nearest_sample=False)
                # Check if sliced trace is empty
                if tr_sliced.data.size == 0:
                    # If so, exit for-loop with break and while-loop with flag
                    flag = False
                    break
                # Trace is not empty -> save it
                print('\nResult of slicing trace with current start time and end time:')
                print(tr_sliced)
                output_filename = (tr_sliced.stats.station + "."
                                   + str(tr_sliced.stats.starttime.year)[-2:]
                                   + ("0" + str(tr_sliced.stats.starttime.month))[-2:]
                                   + ("0" + str(tr_sliced.stats.starttime.day))[-2:]
                                   + ("0" + str(tr_sliced.stats.starttime.hour))[-2:]
                                   + ("0" + str(tr_sliced.stats.starttime.minute))[-2:]
                                   + ("0" + str(tr_sliced.stats.starttime.second))[-2:]
                                   + "." + str(tr_sliced.stats.channel))
                print(f'File name for saving cut trace:\n{output_filename}')
                # Check if a directory for the current output is created
                # Create it if it does not
                cur_output_dir = os.path.join(output_dir, tr_sliced.stats.station)
                if (not os.path.isdir(cur_output_dir)):
                    os.mkdir(cur_output_dir)
                output_path = os.path.join(cur_output_dir, output_filename)
                print(f'Output will be saved to: {output_path}')
                tr_sliced.write(output_path, format="MSEED")

            # Changing the start time / exiting the while-loop
            if cur_starttime.hour < 23:
                cur_starttime.hour += 1
                # Make sure the next piece of data starts
                # with zero minutes and seconds
                cur_starttime.minute = 0
                cur_starttime.second = 0
                cur_starttime.microsecond = 0
            else:
                break
            # cur_starttime = cur_endtime

# # %%
# # flag = True
# st_1 = read('H:/temp/1day/CR10/cr10190808061837.hhz')
# st_2 = read('H:/temp/1day/CR10/cr10190808061942.hhz')
# print(st_1[0].stats)
# print(st_2[0].stats)
