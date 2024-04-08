# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 23:52:53 2021

@author: Igor
"""
import os
import sys
import numpy as np
import pandas as pd
import scipy.signal
from glob import glob
from obspy import read, UTCDateTime
from obspy.core.trace import Trace
from obspy.geodetics import gps2dist_azimuth

def combinations_array(file_list):
    '''
    Creates a two-dimensional array containing in its rows
    the indices of the files from file_list that needs to be correlated.

    Parameters:

    file_list - list of mseed-files;
    the name of each file must be in the STATION.STARTTIME.CHANNEL format

    Returns a 2D numpy array with rows containing the indices of station pair
    to be correlated.
    '''

    # Create a list of tuples that will store all combinations of stations
    combs = []

    # Double for-loop: comparing each file with every other file
    for i in range(len(file_list)):
        # for k in range(1, len(file_list)):
        for k in range(len(file_list)):
            # Extract station name and starttime from file a
            file_name_a = file_list[i].split('\\')[-1]
            station_a = file_name_a.split('.')[0]
            starttime_a = file_name_a.split('.')[1]
            # Extract station name and starttime from file b
            file_name_b = file_list[k].split('\\')[-1]
            station_b = file_name_b.split('.')[0]
            starttime_b = file_name_b.split('.')[1]
            # If the station names do not match
            # and the recording start time is the same,
            # then add the corresponding indexes in combs
            if (station_a != station_b and starttime_a == starttime_b and i < k):
                combs.append((i, k))

    # Before returning, covert list to numpy array
    return np.array(combs)

def combinations_df(file_list, save_excel=''):
    '''
    Creates a pandas DataFrame object containing in its rows
    the indices of the files from file_list that need to be correlated.

    Parameters:

    file_list - list of mseed-files;
    the name of each file must be in the STATION.STARTTIME.CHANNEL format

    save_excel - (optional) path to save the resulting DataFrame object
    in xlsx format

    Returns a complex MultiIndex DataFrame object
    with rows containing the indices of station pair to be correlated.
    '''

    # Create an empty DataFrame that will store all combinations of stations
    df = pd.DataFrame()

    # Double for-loop: comparing each file with every other file
    for i in range(len(file_list)):
        for k in range(len(file_list)):
            # Extract station name and starttime from file a
            file_name_a = file_list[i].split('\\')[-1]
            station_a = file_name_a.split('.')[0]
            # This will extract the entire string corresponding to the starttime
            # starttime_a = file_name_a.split('.')[1]
            # Due to the fact that there may be files that don't start
            # from the first minute of the hour,
            # it makes sense to extract the starttime only up to hours.
            starttime_a = file_name_a.split('.')[1][:-4]

            # Extract station name and starttime from file b
            file_name_b = file_list[k].split('\\')[-1]
            station_b = file_name_b.split('.')[0]
            # This will extract the entire string corresponding to the starttime
            # starttime_b = file_name_b.split('.')[1]
            # Due to the fact that there may be files that don't start
            # from the first minute of the hour,
            # it makes sense to extract the starttime only up to hours.
            starttime_b = file_name_b.split('.')[1][:-4]

            # If the station names do not match
            # and the recording start time is the same
            # and i < k (which allows to exclude repetitions)
            # then add the corresponding indexes to df
            if (station_a != station_b and starttime_a == starttime_b and i < k):
                columns = pd.MultiIndex.from_product(
                        [[f"{station_a}-{station_b}"],
                         [f"{station_a}", f"{station_b}"]])
                utc_time = UTCDateTime('20' + starttime_a)
                index = [str(utc_time)]
                temp_df = pd.DataFrame([[i, k]], columns=columns, index=index)
                df = pd.concat([df, temp_df])

    # Group by indices and sum to eliminate duplicate indices and NaN-values
    df = df.groupby(df.index).sum()

    # Save if necessary
    if save_excel:
        df.to_excel(save_excel)

    return df


def combinations_df_v2(subfolders_files, nsubfolders, subfolders_nfiles,
                       save_excel=''):
    '''
    Creates a pandas DataFrame object containing in its rows
    the indices of the files that need to be correlated.

    Parameters:

    subfolders_files - a list of lists of full paths to files from different
    stations; the name of each file must be in the STATION.STARTTIME.CHANNEL
    format

    nsubfolders - a number of subfolders in the input folder,
    as well as the number of stations used in the calculations

    subfolders_nfiles - a list that conaints the number of files in each
    subfolder

    save_excel - (optional) path to save the resulting DataFrame object
    in xlsx format

    Returns a complex MultiIndex DataFrame object
    with rows containing indices (from variable files) of station pairs
    to be correlated.
    '''

    # Creation of a term for transformation of
    # an index from a separate one for each group of files
    # to a pass-through index for all files
    indx_factor = np.zeros(nsubfolders)
    indx_factor[1:] = np.cumsum(subfolders_nfiles[0:-1])

    # Create an empty DataFrame that will store all combinations of stations
    df = pd.DataFrame()

    # Loop through every folder
    for i in range(nsubfolders):
        # Loop through every file in the current folder
        for j in range(subfolders_nfiles[i]):
        # for j in range(200):
            # Extract station name and starttime from file a
            file_name_a = subfolders_files[i][j].split('\\')[-1]
            station_a = file_name_a.split('.')[0]
            # Due to the fact that there may be files that don't start
            # from the first minute of the hour,
            # it makes sense to extract the starttime only up to hours
            starttime_a = file_name_a.split('.')[1][:-4]
            # Loop through all the subsequent folders from the list
            # This will give all the required pairs of station indices,
            # since if pair a-b is found, pair b-a is not needed
            for k in range(i + 1, nsubfolders):
                # Search the current folder for indices of files
                # with names containing starttime_a
                # NOTE: index_b is tuple; index_b[0] - a numpy array,
                # index_b[0][0] - the first element of a numpy array
                index_b = np.where(
                    np.char.find(subfolders_files[k], starttime_a) > 0)
                # The case if one such file is found
                if index_b[0].size == 1:
                    index_b = index_b[0][0]  # Getting rid of the np.array
                    # Extract station name and starttime (up to hours)
                    # from file b
                    file_name_b = subfolders_files[k][index_b].split('\\')[-1]
                    station_b = file_name_b.split('.')[0]
                    starttime_b = file_name_b.split('.')[1][:-4]
                    print("\nFound pair of stations:")
                    print("Station name\tStarttime")
                    print(f"{station_a}\t{starttime_a}")
                    print(f"{station_b}\t{starttime_b}")
                    # Create a temporary DataFrame for the current files
                    # and concatenate it to the base DataFrame.
                    columns = pd.MultiIndex.from_product(
                            [[f"{station_a}-{station_b}"],
                             [f"{station_a}", f"{station_b}"]])
                    utc_time = UTCDateTime('20' + starttime_a)
                    index = [str(utc_time)]
                    temp_df = pd.DataFrame([[j + indx_factor[i],
                                             index_b + indx_factor[k]]],
                                           columns=columns,
                                           index=index)
                    df = pd.concat([df, temp_df])
                # If more than one pair is found for file_name_a
                elif index_b[0].size > 1:
                    print("\nMore than one pair was found for the current file.")
                    print(f"Current file: {file_name_a}")
                    print("Pairs for the current file:"
                          f"\n{subfolders_files[k][index_b[0]]}")
                    sys.exit()
                # If no pair is found for file A, then the loop continues
                # without doing anything
                else:
                    continue
    # Group by indices and sum to eliminate duplicate indices and NaN-values
    # (Those appear because of how pd.concat operates)
    df = df.groupby(df.index).sum()

    # Save if necessary
    if save_excel:
        df.to_excel(save_excel)

    return df


def combinations_df_v3(files, nsubfolders, subfolders_nfiles,
                       save_excel=''):
    '''
    Creates a pandas DataFrame object containing in its rows
    the indices of the files that need to be correlated.

    Parameters:

    files - an array with full paths to files from different
    stations; the name of each file must be in the STATION.STARTTIME.CHANNEL
    format

    nsubfolders - a number of subfolders in the input folder,
    as well as the number of stations used in the calculations

    subfolders_nfiles - a list that conaints the number of files in each
    subfolder

    save_excel - (optional) path to save the resulting DataFrame object
    in xlsx format

    Returns a complex MultiIndex DataFrame object
    with rows containing indices (from variable files) of station pairs
    to be correlated.
    '''

    # Cumulative sum of files from each subfolder.
    # The elements of this array are the index values of files. They denote
    # the transition from files from one station to files from another station.
    nfiles_cumsum = np.zeros(nsubfolders + 1, dtype='int')
    nfiles_cumsum[1:] = np.cumsum(subfolders_nfiles)

    # Create an empty DataFrame that stores all the file combinations
    # for cross-correlation.
    df = pd.DataFrame()

    # Extract station names from the elements of files array
    stations = np.array([file.split('\\')[-2] for file in files])
    # Extract starttimes from the elements of files array
    # Due to the fact that there may be files that don't start
    # from the first minute of the hour,
    # it makes sense to extract the starttime only up to hours
    times = np.array([file[-16:-8] for file in files])
    # stations[i] and times[i] correspond to one element of files!

    # Loop through files until the last station is reached
    for i in range(nfiles_cumsum[-2]):
    # for i in range(5):
        # Variables that contain various information (for convenience)
        # about the current file (file A)
        # for which pairs are being searched for.
        file_name_a = str(files[i])
        station_a = str(stations[i])
        starttime_a = str(times[i])
        print(f"\nThe current file: \n{file_name_a}")

        # This line determines the index of files,
        # from which the search for pairs for file A will be performed.
        # This is done by finding the index of the first value of nfiles_cumsum
        # that is greater than i (index of file A).
        index_factor = np.argmax(i < nfiles_cumsum)
        # The next line searches for files B from other stations
        # with the same starttime as file A.
        # We shift the indexes of the times array by the index_factor value
        # so that pairs to file A are only sought from other stations
        # further down in the list.
        index_b = np.where(
            np.char.find(times[nfiles_cumsum[index_factor]:], starttime_a) == 0)

        # Now there are only two possibilities:
        # either pairs for file A were found and index_b[0] size is not equal
        # to zero,
        # or pairs were not found and index_b[0] size is equal to zero.

        # Case 1: some pairs are found
        if index_b[0].size != 0:
            # Arrays containing information about files B -
            # pairs of file A(for convenience)
            file_names_b = files[index_b[0][:] + nfiles_cumsum[index_factor]]
            stations_b = stations[index_b[0][:] + nfiles_cumsum[index_factor]]
            # this line will be relevant later
            print("Pairs to the current file:")

            # There is a possibility that files B will contain several files
            # from the same station. This is an error.

            # Case 1: there are no repeating stations in the files B
            if len(np.unique(stations_b)) == len(stations_b):
                # Create a temporary DataFrame for the current files A and b;
                # 1. Create a pd.MultiIndex for columns of a temporary DataFrame
                list_of_tuples = []
                for j in range(len(file_names_b)):
                    # Print the name of the current files B
                    print(f"{file_names_b[j]}")
                    level_0 = f"{station_a}-{stations_b[j]}"
                    level_1_0 = f"{station_a}"
                    level_1_1 = f"{stations_b[j]}"
                    list_of_tuples.append((level_0, level_1_0))
                    list_of_tuples.append((level_0, level_1_1))
                columns = pd.MultiIndex.from_tuples(list_of_tuples)
                # 2. Create index of this temporary DataFrame
                utc_time = UTCDateTime('20' + starttime_a)
                index = [str(utc_time)]
                # 3. Create the data vector itself
                data = np.full(shape=(1, 2*len(index_b[0])), fill_value=i,
                               dtype=int)
                data[0, 1::2] = index_b[0] + nfiles_cumsum[index_factor]
                # 4. Put everything together in a DataFrame object.
                temp_df = pd.DataFrame(data,
                                       columns=columns,
                                       index=index)
                # 5. Concatenation of a temporary DataFrame
                # to a "total" DataFrame
                df = pd.concat([df, temp_df])
            # Case 2: some of the files B contain repeating stations
            else:
                # Print names of the files B and display an error.
                for j in range(len(file_names_b)):
                    print(f"{file_names_b[j]}")
                print('Multiple files from the same station are detected.' +
                      ' Error!')
                sys.exit()
        # Case 1: no pairs are found
        else:
            # If no pairs are found for file A, then the loop continues
            # without doing anything
            print("There are no pairs to the current file!")
            continue

    # Group by indices and sum to eliminate duplicate indices and NaN-values
    # (Those appear because of how pd.concat operates)
    df = df.groupby(df.index).sum()

    # Save if necessary
    if save_excel:
        df.to_excel(save_excel)

    return df


def corr_length_samples(length_sec, file_list):
    '''
    Calculates the required cross-correlation length in sample counts
    based on the transmitted cross-correlation length in seconds.

    Parameters:

    length_sec - cross-correlation length in seconds

    file_list - list of mseed-files;
    the name of each file must be in the STATION.STARTTIME.CHANNEL format

    Returns cross-correlation length in sample counts.
    '''
    # Define variables
    st = read(file_list[0])
    srate = int(st[0].stats.sampling_rate)

    # Calculate the length of cross-correlation in samples
    # length_samples = int(2*length_sec*srate + 1)
    length_samples = int(length_sec*srate)

    return length_samples


def corr_zerolag_equal(file_list):
    '''
    Finds the index corresponding to the zero cross-correlation lag
    (When the centers of the correlated signals are aligned).
    Important note: we assume that all data are the same length.

    Parameters:

    file_list - list of mseed-files;
    the name of each file must be in the STATION.STARTTIME.CHANNEL format

    Return the index of the zero lag.
    '''
    st = read(file_list[0])
    lags = scipy.signal.correlation_lags(len(st[0].data), len(st[0].data),
                                         mode='same')
    zero_lag = np.argmin(lags**2)

    return zero_lag


def corr_zerolag(sig1, sig2):
    '''
    Finds the index corresponding to the zero cross-correlation lag
    (When the centers of the correlated signals are aligned).

    Parameters:

    sig1, sig2 - Numpy arrays containing data from two different stations.
    They can be of different lengths.

    Return the index of the zero lag.
    '''

    # Reminder - 'mode' argument:
    # mode='full' - full discrete linear cross-correlation of the inputs,
    # lenght n - k - 1
    # mode='valid' - the output consists only of those elements that
    # do not rely on the zero-padding
    # mode='same' - the output is the same size as in1,
    # centered with respect to the ‘full’ output.

    lags = scipy.signal.correlation_lags(len(sig1), len(sig2), mode='full')
    # Zero lag - when zero value of sig 1 is alligned with zero value of sig2
    zero_lag = np.argmin(lags**2)

    return zero_lag

def corr_days(input_folder):
    '''
    Analyzes folders containing data from individual stations
    and determines the duration (in days) of the longest dataset.

    Parameters:

    input_folder - root working directory containing folders
    with one hour long mseed files from various stations

    Returns the number of days recorded
    by the station that was operational the longest.
    '''

    # Determine how many folders (stations) the input_folder contains.
    folders = glob(input_folder + '/*')
    num_files = np.zeros(len(folders))

    # Loop through each folder and count the files inside.
    for i in range(len(folders)):
        num_files[i] = len(os.listdir(folders[i]))
    # Determine the data for how many days are contained in the folder
    # with the most files.
    days = int(np.ceil(np.max(num_files)/24))

    return days

def load_signals(file_list, i, j, cons_len=False):
    '''
    Loads data from MSEED files containing only one trace and saves
    stream objects and signal data arrays into variables sorted according to
    load order (i - first, j - second) or signal length.

    Parameters:

    file_list - list containing strings describing the full paths
    to the processed files

    i, j - file_list indices corresponding to a pair of files

    cons_len - whether to sort the result depending on the length of the signal

    Returns two pairs of sorted Obspy streams and Numpy arrays.
    '''
    # Load the first and second streams
    temp_st1 = read(file_list[i])
    temp_st2 = read(file_list[j])

    if cons_len:
        # If the length of the two datasets is different,
        # determine which one is longer
        # The longer one will be assigned to sig1,
        # the shorter one will be assigned to sig2
        # If the length is the same, then st1 = sig1, st2 = sig2
        if len(temp_st1[0].data) >= len(temp_st2[0].data):
            # 1
            st1 = temp_st1
            sig1 = st1[0].data
            # 2
            st2 = temp_st2
            sig2 = st2[0].data
        elif len(temp_st1[0].data) < len(temp_st2[0].data):
            # 1
            st1 = temp_st2
            sig1 = st1[0].data
            # 2
            st2 = temp_st1
            sig2 = st2[0].data
        # To be on the safe side, return both numbered signals and stream objects
        return st1, sig1, st2, sig2

    return temp_st1, temp_st1[0].data, temp_st2, temp_st2[0].data


def corr_data_to_trace(data, station1, station2, stations_info, stats_dict,
                       sac=False):
    '''
    Creates a Trace instance of the Obspy package from
    the calculated cross-correlations and various metadata.

    Parameters:

    data - numpy array containing the calculated cross-correlation

    station1, station2 - strings containing the names of two stations
    for which the cross-correlation was calculated

    stations_info - DataFrame instance containing information about
    the latitude, longitude and elevation of all stations involved
    in processing

    stats_dict - Python dictionary, the keys of which are
    the names of stations, and the values contain metadata obtained
    from one of the processed files from the corresponding station

    sac - whether to add the metadata required for correct
    saving in the SAC format.

    Returns a Trace object of the Obspy package with correctly filled headers
    and calculated cross-correlation.
    '''
    # Create Obspy trace object with data but no metadata (yet)
    tr = Trace(data)

    # Filling in the necessary metadata
    tr.stats['network'] = (f"{stats_dict[station1].network}"
                           f"-{stats_dict[station2].network}")
    tr.stats['station'] = station1 + station2
    tr.stats['location'] = '-'
    tr.stats['channel'] = (f"{stats_dict[station1].channel}"
                           f"-{stats_dict[station2].channel}")
    tr.stats['sampling_rate'] = stats_dict[station1].sampling_rate

    # If output will be eventually saved in SAC format, these fields must be
    # filled as well
    if sac:
        tr.stats['sac'] = {}
        tr.stats['sac']['stla'] = stations_info.loc[station1]['latitude']
        tr.stats['sac']['stlo'] = stations_info.loc[station1]['longitude']
        tr.stats['sac']['stel'] = stations_info.loc[station1]['elevation']
        tr.stats['sac']['evla'] = stations_info.loc[station2]['latitude']
        tr.stats['sac']['evlo'] = stations_info.loc[station2]['longitude']
        tr.stats['sac']['evel'] = stations_info.loc[station2]['elevation']

        # get distance, azimuth and back azimuth between stations
        dist, az, baz = gps2dist_azimuth(tr.stats.sac.stla, tr.stats.sac.stlo,
                                         tr.stats.sac.evla, tr.stats.sac.evlo)
        tr.stats['sac']['dist'] = dist / 1000
        tr.stats['sac']['az'] = az
        tr.stats['sac']['baz'] = baz

    return tr


def create_stats_dict(input_folder):
    '''
    Creates a Python dictionary, the keys of which are the names of the stations
    and the values are the corresponding Stats objects (obspy module).
    Note: file structure required is input_folder/station_name/data;
    it is assumed that the metadata of the stations does not change

    Parameters:

    input_folder - a string with the full path to the working folder

    Returns a Python dictionary containing station metadata
    '''

    # Define some variables
    stats_dict = {}
    st = read()
    folder_list = glob(input_folder + '/*/')

    # Loop through every folder in folder_list...
    for folder in folder_list:
        # ... and load first file in it
        print(os.listdir(folder)[0])
        st = read(folder + os.listdir(folder)[0])
        # Save its name (key) and stats (value) in the dictionary
        stats_dict[st[0].stats.station] = st[0].stats

    return stats_dict
