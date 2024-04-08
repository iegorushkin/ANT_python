import shutil, os
import glob
import numpy as np
import matplotlib.pyplot as plt
from obspy import read, UTCDateTime


def check_for_duplicates(x):
    '''
    Searches for duplicate values in numpy array x
    and returns a Boolean array labeled True in place of duplicates.
    '''
    idup = np.zeros(len(x), dtype='bool_')

    for i in range(len(x) - 1):
        if x[i] == x[i + 1]:
            idup[i], idup[i + 1] = True, True
        else:
            continue

    return idup


# Setup

input_dir = r'H:\Central_Kamchatka\Raw_data\Permanent_stations'
output_dir = r'H:\CK_RAW_Z\1day\perm_stations'

station_list_temp = os.listdir(input_dir)
station_list = []
for station in station_list_temp:
    if station[-2:] == 'BH':
        station_list.append(station)
del station_list_temp

# Calculations

for station in station_list:
    # infilist1 = os.listdir(input_dir)

    # Check if the output directory for the current station is created
    # Create it if it does not
    cur_output_dir = os.path.join(output_dir, station)
    if (not os.path.isdir(cur_output_dir)):
        os.mkdir(cur_output_dir)

    # File paths corresponding to the current station and the desired channel
    file_list = np.array(glob.glob(input_dir + '/' + station + '/*z'))

    # Extraction of dates from the names of each file
    dates = np.empty(len(file_list), dtype='U10')
    for i in range(len(file_list)):
        dates[i] = file_list[i][-16:-10]
        # dates[i] = file_list[i][-20:-10]

    # Find the indices of the dates that repeat and the number of those repeats
    # If the date repeats, then the corresponding files need to be merged.
    i_torn_days = check_for_duplicates(dates)
    _, counts = np.unique(dates[i_torn_days], return_counts=True)

    # Paths to files that just need to be copied
    files2move = file_list[~i_torn_days]
    # Paths to files that need to be merged
    files2merge = file_list[i_torn_days]
    # Thanks to this offset, the loop will iterate over all file paths
    # from the variable files2merge
    ind_offset = 0

    for i, j in enumerate(counts):
        print(f'\nGroup number {i + 1}')
        # Create empty stream
        st = read()
        st.clear()

        # Each loop - files corresponding to a specific day
        for k in range(j):
            print(f"{files2merge[k + ind_offset]}")
            st += read(files2merge[k + ind_offset])

        # All overlapped samples and gaps are linearly interpolated
        st.merge(method=1, interpolation_samples=-1, fill_value='interpolate')
        # Save this merged stream
        output_filename = files2merge[ind_offset][-20:]
        output_path = os.path.join(cur_output_dir, output_filename)
        print(f'Output will be saved to: {output_path}')
        st.write(output_path, format="MSEED")

        ind_offset += j

    # Copy files that don't need to be merged to a new folder.
    for file in files2move:
        print(f"\nCopying {file}...")
        shutil.copy(file, cur_output_dir)
        print("Done.")
# %%

# import matplotlib.pyplot as plt

# st1 = read(files2merge[0])
# st2 = read(files2merge[1])
# st = (st1 + st2)

# # sort
# st.sort(['starttime'])
# # start time in plot equals 0
# dt = st[0].stats.starttime.timestamp

# fig, ax = plt.subplots(3, 1, figsize=(10, 6), sharex=True)

# for i in range(2):
#     t = np.linspace(st[i].stats.starttime.timestamp - dt,
#                     st[i].stats.endtime.timestamp - dt,
#                     st[i].stats.npts)
#     ax[i].plot(t, st[i].data)

# st.merge(method=1, interpolation_samples=-1, fill_value=0)
# t = np.linspace(st[0].stats.starttime.timestamp - dt,
#                 st[0].stats.endtime.timestamp - dt,
#                 st[0].stats.npts)
# ax[2].plot(t, st[0].data)
