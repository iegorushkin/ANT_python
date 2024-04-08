import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.fftpack
import os
import sys
from scipy import signal, stats


def visualize_crosscorr(cc_matrix, file_list, dont_open=True, save_path=''):
    '''
    Performs neat visualization of the correlation matrix
    and averaged cross-correlation for a pair of stations.

    Parameters:

    cc_matrix - pandas DataFrame object, each column of which is corresponding
    to the mean cross-correlation of all data from a pair of stations

    file_list - list of strings describing the entire path to xlsx files
    that contain daily cross correlations for different pairs of stations

    dont_open - display figures or not

    save_path - path to the folder where the pictures will be saved
    (if required to do so)
    '''

    for i in range(len(file_list)):
        # Load current xlsx
        cc_matrix_days = pd.read_excel(file_list[i], index_col=0,
                                       header=(0, 1))
        # Plot!
        fig, ax = plt.subplots(2, 1, figsize=(22, 10),  sharex=True,
                               gridspec_kw={'height_ratios': [2, 1]},
                               constrained_layout=True)
        if dont_open:
            plt.close(fig)  # So that the figure is not shown
        fig.suptitle("Cross-correlations of records from stations" +
                     f" {cc_matrix.columns[i]}", fontsize=16)

        # Subplot 1
        cb = ax[0].contourf(cc_matrix_days.index,
                            cc_matrix_days.columns.levels[1],
                            cc_matrix_days.T,
                            cmap=plt.cm.turbo, levels=21)
        fig.colorbar(cb, ax=ax[0], orientation='vertical',
                     label='Amplitude, a.u.', shrink=0.75, aspect=20, pad=0.01)
        ax[0].set_ylabel('Day', fontsize=12)
        # Deal with y-ticks:
        # Make them more sparse and rotate the labels
        y_values = np.arange(len(cc_matrix_days.columns.levels[1]))[::20]
        y_ticks = cc_matrix_days.columns.levels[1][::20]
        ax[0].set_yticks(y_values)
        ax[0].set_yticklabels(y_ticks, rotation=45)
        # TO DO: make the number of labels on the y-axis more standardized.

        # Subplot 2
        ax[1].plot(cc_matrix.index, cc_matrix.iloc[:, i], linewidth=2)
        ax[1].set_title('Average of the daily cross-correlation functions',
                        fontsize=14)
        ax[1].set_xlabel('Lag, sec', fontsize=12)
        ax[1].set_ylabel('Amplitude, a.u.', fontsize=12)

        # Save figure if required
        if save_path:
            file_name = cc_matrix.columns[i]
            flag = 1
            while True:
                if os.path.exists(save_path+file_name+'_'+str(flag)+'.png'):
                    flag += 1
                else:
                    break
            fig_file_name = save_path + file_name + '_' + str(flag) + '.png'
            fig.savefig(fig_file_name, bbox_inches='tight', dpi=100)
            print("Saved to " + fig_file_name)


def visualize_signal(st, freq_xlim=[0, 1], dont_open=True, file_name='',
                     save_path='', time_xlim=[0, 60]):
    '''
    Plots a single trace from the obspy stream object st
    in the time and frequency domain.

    Parameters:

    st - obspy stream object containing only 1 trace

    file_name - the name of the mseed file from which the stream is read

    freq_xlim - range of frequencies to plot

    save_path - path to the folder where the pictures will be saved
    (if required to do so)

    dont_open - display figures or not
    '''

    # Define variables for calculations and plotting
    sig = st[0].data
    npts = len(sig)
    srate = int(st[0].stats.sampling_rate)
    time = np.arange(npts)/(srate*60)  # Time in minutes
    hz = np.linspace(0, srate/2, int(np.floor(npts/2 + 1)))

    # Amplitude spectrum of the original signal
    amp_sig = np.abs(scipy.fftpack.fft(sig)/npts)
    amp_sig = amp_sig[:len(hz)]
    amp_sig[1:] = 2*amp_sig[1:]

    # Create fig and ax objects
    fig1, ax1 = plt.subplots(2, 1, figsize=(25.0, 50.0), dpi=110)
    # fig1.tight_layout()
    if dont_open:
        plt.close(fig1)  # So that the figure is not shown

    # Title of the figure
    fig1.suptitle(f"Station: {st[0].stats.station}, "
                  f"channel: {st[0].stats.channel}."
                  f"\nDate: {st[0].stats.starttime.day}/"
                  f"{st[0].stats.starttime.month}/"
                  f"{st[0].stats.starttime.year}."
                  f" Time: {str(st[0].stats.starttime)[-16:-8]}"
                  f" - {str(st[0].stats.endtime)[-16:-8]}",
                  fontsize=16, y=0.96)

    # time domain
    ax1[0].plot(time, sig, 'k', label='signal')
    ax1[0].set_xlabel('Time (min)', fontsize=12)
    ax1[0].set_ylabel('Amplitude', fontsize=12)
    ax1[0].set_title("Time domain", fontsize=14)
    ax1[0].set_ylim([np.min(sig), np.max(sig)])
    ax1[0].set_xlim(time_xlim)
    # ax1[0].legend(fontsize=12)
    ax1[0].grid()
    # frequency domain
    ax1[1].plot(hz, amp_sig, 'k', label='signal')
    ax1[1].set_xlabel('Frequency (Hz)', fontsize=12)
    ax1[1].set_ylabel('Amplitude', fontsize=12)
    ax1[1].set_title("Frequency domain", fontsize=14)
    ax1[1].set_ylim([0, 1.1*np.max(amp_sig[1:])])
    ax1[1].set_xlim(freq_xlim)
    # ax1[1].legend(fontsize=12)
    ax1[1].grid()

    # Save figure if required
    if save_path and file_name:
        flag = 1
        while True:
            if os.path.exists(save_path+file_name+'_'+str(flag)+'.png'):
                flag += 1
            else:
                break
        fig_file_name = save_path + file_name + '_' + str(flag) + '.png'
        fig1.savefig(fig_file_name, bbox_inches='tight')
        print("Saved to " + fig_file_name)

    elif save_path or file_name:
        print("Both save_path and file_name are required to save the picture.")


def visualize_signal_adv(st, freq_xlim=[0, 1], dont_open=True, file_name='',
                         save_path='',):
    '''
    Plots a single trace from the obspy stream object st
    in the time and frequency domain. Removes mean and linear trend,
    then also plots this result.

    Parameters:

    st - obspy stream object containing only 1 trace

    file_name - the name of the mseed file from which the stream is read

    freq_xlim - range of frequencies to plot

    save_path - path to the folder where the pictures will be saved
    (if required to do so)

    dont_open - display figures or not
    '''

    # save_path = 'E:/Work/Projects/Central_Kamchatka/Pictures/'

    # Define variables for calculations and plotting
    sig = st[0].data
    npts = len(sig)
    srate = int(st[0].stats.sampling_rate)
    time = np.arange(npts)/(srate*60)  # Time in minutes
    hz = np.linspace(0, srate/2, int(np.floor(npts/2 + 1)))

    # Amplitude spectrum of the original signal
    amp_sig = np.abs(scipy.fftpack.fft(sig)/npts)
    amp_sig = amp_sig[:len(hz)]
    amp_sig[1:] = 2*amp_sig[1:]

    # Create fig and ax objects
    fig1, ax1 = plt.subplots(2, 2, figsize=(25.0, 50.0), dpi=110)
    # fig1.tight_layout()
    if dont_open:
        plt.close(fig1)  # So that the figure is not shown

    # Title of the figure
    fig1.suptitle(f"Station: {st[0].stats.station}, "
                  f"channel: {st[0].stats.channel}."
                  f"\nDate: {st[0].stats.starttime.day}/"
                  f"{st[0].stats.starttime.month}/"
                  f"{st[0].stats.starttime.year}."
                  f" Time: {str(st[0].stats.starttime)[-16:-8]}"
                  f" - {str(st[0].stats.endtime)[-16:-8]}",
                  fontsize=16, y=0.96)

    # Plot the unchanged signal
    # time domain
    ax1[0, 0].plot(time, sig, 'k', label='signal')
    ax1[0, 0].set_xlabel('Time (min)', fontsize=12)
    ax1[0, 0].set_ylabel('Amplitude', fontsize=12)
    ax1[0, 0].set_title("Time domain", fontsize=14)
    ax1[0, 0].set_xlim([0, 60])
    ax1[0, 0].legend(fontsize=12)
    ax1[0, 0].grid()
    # frequency domain
    ax1[0, 1].plot(hz, amp_sig, 'k', label='amp. spectrum of the signal')
    ax1[0, 1].set_xlabel('Frequency (Hz)', fontsize=12)
    ax1[0, 1].set_ylabel('Amplitude', fontsize=12)
    ax1[0, 1].set_title("Frequency domain", fontsize=14)
    ax1[0, 1].set_xlim(freq_xlim)
    ax1[0, 1].legend(fontsize=10)
    ax1[0, 1].grid()

    # Remove mean and linear trend from the signal
    msig = sig - np.mean(sig)
    scipy.signal.detrend(msig, type='linear', overwrite_data=True)

    # Amplitude spectrum
    amp_msig = np.abs(scipy.fftpack.fft(msig)/npts)
    amp_msig = amp_msig[:len(hz)]
    amp_msig[1:] = 2*amp_msig[1:]

    # plot the result.
    # time domain
    ax1[1, 0].plot(time, msig, 'r',
                   label='signal after removing\nmean and trend')
    ax1[1, 0].set_xlabel('Time (min)', fontsize=12)
    ax1[1, 0].set_ylabel('Amplitude', fontsize=12)
    # ax1[1, 0].set_title("Time domain")
    ax1[1, 0].set_xlim([0, 60])
    ax1[1, 0].legend(fontsize=12)
    ax1[1, 0].grid()
    # frequency domain
    ax1[1, 1].plot(hz, amp_msig, 'r',
                   label='amp. spectrum of the signal\nafter removing mean and trend')
    ax1[1, 1].set_xlabel('Frequency (Hz)', fontsize=12)
    ax1[1, 1].set_ylabel('Amplitude', fontsize=12)
    # ax1[1, 1].set_title("Frequency domain")
    ax1[1, 1].set_xlim(freq_xlim)
    ax1[1, 1].legend(fontsize=10)
    ax1[1, 1].grid()

    # Setting up the y-axis for time-domain plots
    ymin_t = np.min([np.min(sig), np.min(msig)])
    ymax_t = np.max([np.max(sig), np.max(msig)])
    ax1[0, 0].set_ylim([ymin_t, ymax_t])
    ax1[1, 0].set_ylim([ymin_t, ymax_t])

    # Setting up the y-axis for frequency-domain plots
    ymin_f = 0
    ymax_f = np.min([np.max(amp_sig[:len(hz)]),
                     np.max(amp_msig[:len(hz)])])
    ax1[0, 1].set_ylim([ymin_f, ymax_f])
    ax1[1, 1].set_ylim([ymin_f, ymax_f])

    # Save figure if required
    if save_path and file_name:
        flag = 1
        while True:
            if os.path.exists(save_path+file_name+'_'+str(flag)+'.png'):
                flag += 1
            else:
                break
        fig_file_name = save_path + file_name + '_' + str(flag) + '.png'
        fig1.savefig(fig_file_name, bbox_inches='tight')
        print("Saved to " + fig_file_name)

    elif save_path or file_name:
        print("Both save_path and file_name are required to save the picture.")


def visualize_fsignal(st, frange, hp_order, lp_order, freq_xlim=[0, 2],
                      dont_open=True, file_name='', save_path=''):
    '''
    Plots the normalized bandpass filtered signal (contained in first trace
    of the obspy object st) and the frequency response
    of the used high- and lowpass filters.

    Parameters:

    st - obspy stream object that contains 1 filtered trace

    frange - bandwidth of the combined high- and lowpass filter
    used on the signal

    hp_order - the order of the highpass filter used on the signal

    lp_order - the order of the lowpass filter used on the signal

    dont_open - display figures or not

    file_name - the name of the mseed file from which the stream is read

    save_path - path to the folder where the pictures will be saved
    (if required to do so)
    '''
    # Define variables for calculations, plotting and convenience
    sig = st[0].data
    npts = len(sig)
    srate = int(st[0].stats.sampling_rate)
    hz = np.linspace(0, srate/2, int(np.floor(npts/2 + 1)))

    # Ideal shape of the combined filter
    filter_xshape = [0, frange[0], frange[0], frange[1], frange[1], srate/2]
    filter_yshape = [0, 0, 1, 1, 0, 0]

    # generate the impulse
    impulse = np.zeros(npts)
    impulse[int(npts/2)] = 1

    # create coefficients for 2 Butter filters - highpass and lowpass
    hp_b, hp_a = signal.butter(hp_order, frange[0], btype='highpass', fs=srate)
    lp_b, lp_a = signal.butter(lp_order, frange[1], btype='lowpass', fs=srate)

    # apply these filters to the impulse.
    hp_impulse = signal.lfilter(hp_b, hp_a, impulse)
    lp_impulse = signal.lfilter(lp_b, lp_a, impulse)

    # compute amplitude spectra of filtered impulses
    hp_impulse_amp = np.abs(scipy.fftpack.fft(hp_impulse))
    hp_impulse_amp = hp_impulse_amp[0:len(hz)]
    lp_impulse_amp = np.abs(scipy.fftpack.fft(lp_impulse))
    lp_impulse_amp = lp_impulse_amp[0:len(hz)]

    combined_impulse_amp = np.zeros(len(hz))
    # combine spectra of the two filters into one
    # 1) find the mean value of the frequency range
    mean_freq = np.mean(frange)
    # 2) find which index this value corresponds to
    ind = np.argmin((hz - mean_freq)**2)
    # 3) The values for the combined spectrum before this index
    # are taken from the spectrum of the highpass filter, the values after
    # this index are taken from the spectrum of the lowpass filter
    combined_impulse_amp[:ind] = hp_impulse_amp[:ind]
    combined_impulse_amp[ind:] = lp_impulse_amp[ind:]

    # Amplitude spectrum of the filtered signal
    amp_sig = np.abs(scipy.fftpack.fft(sig)/npts)
    amp_sig = amp_sig[:len(hz)]
    amp_sig[1:] = 2*amp_sig[1:]

    # Spectrum normalization
    amp_sig = amp_sig / np.max(amp_sig)

    # Plot: 1) spectrum of the filtered signal, 2) ideal filter shape,
    # 3) real filter shape
    fig, ax = plt.subplots(figsize=(25.0, 50.0), dpi=110)

    if dont_open:
        plt.close(fig)  # So that the figure is not shown

    # Title of the figure
    fig.suptitle(f"Station: {st[0].stats.station}, "
                 f"channel: {st[0].stats.channel}."
                 f"\nDate: {st[0].stats.starttime.day}/"
                 f"{st[0].stats.starttime.month}/"
                 f"{st[0].stats.starttime.year}."
                 f" Time: {str(st[0].stats.starttime)[-16:-8]}"
                 f" - {str(st[0].stats.endtime)[-16:-8]}",
                 fontsize=16, y=0.96)

    ax.plot(hz, amp_sig, 'k', label='filtered signal')
    ax.plot(filter_xshape, filter_yshape, 'r', label='ideal filter shape')
    ax.plot(hz, combined_impulse_amp, 'm', label='real filter shape')
    ax.set_xlabel('Frequency (Hz)', fontsize=12)
    ax.set_ylabel('Normalized amplitude', fontsize=12)
    ax.set_title("Frequency domain", fontsize=14)
    ax.set_ylim([0, 1.05])
    ax.set_xlim(freq_xlim)
    ax.legend(fontsize=12)
    ax.grid()

    # Save figure if required
    if save_path and file_name:
        flag = 1
        while True:
            if os.path.exists(save_path+file_name+'_'+str(flag)+'.png'):
                flag += 1
            else:
                break
        fig_file_name = save_path + file_name + '_' + str(flag) + '.png'
        fig.savefig(fig_file_name, bbox_inches='tight')
        print("Saved to " + fig_file_name)

    elif save_path or file_name:
        print("Both save_path and file_name are required to save the picture.")


def visualize_iir_filter(b_coeffs, a_coeffs,
                         filter_name, filter_order, filter_type,
                         filter_xshape, filter_yshape,
                         npts, hz, freq_xlim=[0, 0.1]):
    '''
    Plots the frequency response of an IIR-filter in the frequency - gain and
    frequency - dB coordinates.

    Parameters:

    b_coeffs, a_coeffs - filter kernel coefficients

    filter_name - name of the filter, for example "highpass Butterworth filter"

    filter_order - the order of the filter

    filter_type - 0 or 1; 0 - highpass filter, 1 - lowpass filter

    filter_xshape, filter_yshape - ideal shape of the filter in
    the frequency domain (x - frequency, y - attenuation)

    npts - number of points to generate impulse
    hz - frequency vector that are used with this impulse

    freq_xlim - limit of the x-axis (frequency) on the frequency-domain plots
    '''

    # generate the impulse
    impulse = np.zeros(npts)
    impulse[int(npts/2)] = 1

    # apply the filter
    # fimpulse = signal.filtfilt(b_coeffs, a_coeffs, impulse)
    fimpulse = signal.lfilter(b_coeffs, a_coeffs, impulse)

    # compute amplitude spectrum
    fimpulse_amp = np.abs(scipy.fftpack.fft(fimpulse))
    fimpulse_amp = fimpulse_amp[0:len(hz)]

    # compute power spectrum
    fimpulse_pow = fimpulse_amp**2
    # deal with zero before np.log
    zeroes = np.nonzero(fimpulse_pow == 0)[0]
    for i in zeroes:
        try:
            fimpulse_pow[i] = 0.5*(fimpulse_pow[i - 1]
                                   + fimpulse_pow[i + 1])
        except IndexError:
            fimpulse_pow[i] = 0.5*(fimpulse_pow[i - 1]
                                   + fimpulse_pow[i - 2])

    # compute apltidute spectrum in dB
    fimpulse_db = 10*np.log10(fimpulse_pow)

    # plot the frequency response of the filter
    fig1, ax1 = plt.subplots(2, 1, figsize=(25, 50))
    fig1.suptitle('Frequency response of the ' + filter_name + ', order: '
                  + str(filter_order),
                  fontsize=16,
                  y=0.92)
    # freq. domain (linear)
    ax1[0].plot(filter_xshape, filter_yshape, 'k', label='Ideal', linewidth=2)
    ax1[0].plot(hz, fimpulse_amp, 'r', label='Actual', linewidth=3)
    ax1[0].set_xlim(freq_xlim)
    ax1[0].set_xlabel('Frequency (Hz)', fontsize=12)
    ax1[0].set_ylabel('Gain (percent)', fontsize=12)
    ax1[0].grid()
    ax1[0].legend(fontsize=12)
    # freq. domain (log)
    if filter_type == 0:
        ax1[1].plot(filter_xshape,
                    [np.min(fimpulse_db), np.min(fimpulse_db), 0, 0],
                    'k', label='Ideal', linewidth=2)
    elif filter_type == 1:
        ax1[1].plot(filter_xshape,
                    [0, 0, np.nanmin(fimpulse_db), np.nanmin(fimpulse_db)],
                    'k', label='Ideal', linewidth=2)
        ax1[1].set_ylim([-80, 8])
    ax1[1].plot(hz, fimpulse_db, 'r', label='Actual', linewidth=3)
    ax1[1].set_xlim(freq_xlim)
    ax1[1].set_xlabel('Frequency (Hz)', fontsize=12)
    ax1[1].set_ylabel('Gain (db)', fontsize=12)
    ax1[1].legend(fontsize=12)
    ax1[1].grid()

    # fig1.tight_layout()


def visualize_fir_filter(filter_kernel, filter_name, filter_order, filter_type,
                         filter_xshape, filter_yshape, srate,
                         freq_xlim=[0, 0.1]):
    '''
    Plots the frequency response of an FIR-filter in the frequency - gain and
    frequency - dB coordinates.

    Parameters:

    filter_kernel - filter kernel coefficients

    filter_name - name of the filter, for example "highpass Butterworth filter"

    filter_order - the order of the filter

    filter_type - 0 or 1; 0 - highpass filter, 1 - lowpass filter

    filter_xshape, filter_yshape - ideal shape of the filter in
    the frequency domain (x - frequency, y - attenuation)

    srate - sampling rate of the system

    freq_xlim - limit of the x-axis (frequency) on the frequency-domain plot
    '''
    # compute the frequencies vector for this filter kernel
    filter_hz = np.linspace(0, srate/2, int(np.floor(len(filter_kernel)/2)+1))

    # compute the amplitude spectrum of the filter kernel
    filter_amp = np.abs(scipy.fftpack.fft(filter_kernel))
    filter_amp = filter_amp[0:len(filter_hz)]

    # compute the power spectrum of the filter kernel
    filter_pow = filter_amp**2

    # compute apltidute spectrum in dB
    filter_db = 10*np.log10(filter_pow)

    # plot "kernel" in the frequency domain
    fig1, ax1 = plt.subplots(2, 1, figsize=(25, 50))
    fig1.suptitle('Frequency response of the ' + filter_name + ', order: '
                  + str(filter_order),
                  fontsize=16,
                  y=0.92)
    # freq. domain (linear)
    ax1[0].plot(filter_xshape, filter_yshape, 'k', label='Ideal', linewidth=2)
    ax1[0].plot(filter_hz, filter_amp, 'r', label='Actual', linewidth=3)
    ax1[0].set_xlim(freq_xlim)
    ax1[0].set_xlabel('Frequency (Hz)', fontsize=12)
    ax1[0].set_ylabel('Gain (percent)', fontsize=12)
    ax1[0].legend(fontsize=12)
    ax1[0].grid()
    # freq. domain (log)
    if filter_type == 0:
        ax1[1].plot(filter_xshape, [np.min(filter_db), np.min(filter_db), 0, 0],
                    'k', label='Ideal', linewidth=2)
    elif filter_type == 1:
        ax1[1].plot(filter_xshape,
                    [0, 0, np.min(filter_db), np.min(filter_db)],
                    'k', label='Ideal', linewidth=2)
        ax1[1].set_ylim([-80, 8])
    ax1[1].plot(filter_hz, filter_db, 'r', label='Actual', linewidth=3)
    ax1[1].set_xlim(freq_xlim)
    ax1[1].set_xlabel('Frequency (Hz)', fontsize=12)
    ax1[1].set_ylabel('Gain (db)', fontsize=12)
    ax1[1].legend(fontsize=12)
    ax1[1].grid()

    # fig1.tight_layout()


def visualize_velocity_deviations(coastlines_list, x_limits, y_limits,
                                  stations_path, rays_path,
                                  step=1, yx_aspect_ratio=1, key=0,
                                  dont_open=0, save_path=0):
    '''
    For a given set of periods, creates maps with coastlines, seismic stations
    and rays between them.
    The rays are colored according to the Rayleigh group velocity.
    NOTE: the behavior when one subplot is placed on a Figure is unscripted.

    Parameters:

    coastlines_list - A list of coastlines
    that lie within a previously specified range of geographic coordinates

    x_limits, y_limits - the geographic coordinates within which maps are plotted

    stations_path - the path to the excel-file with information
    about the relevant stations

    rays_path - the path to the directory containing the ray files

    step - visualize every step-file in rays_path

    yx_aspect_ratio - the ratio of the length of the unit segment of the y-axis
    to the unit segment of the x-axis;
    by default, these segments are equal in length

    key - 0 for Figures with 2x2 subplots, 1 for Figures with 3x3 subplots

    dont_open - do not open or open windows with built visualizations;
    by default, the windows will open

    save_path - the path to the directory where Figures will be saved;
    by default, Figures are not saved
    '''

    # Create DataFrame with information about the relevant stations
    stations = pd.read_excel(stations_path, index_col=0)

    # Check if step is legit
    if type(step) != int:
        sys.exit("step argument must be an integer!")
    # Create a list with the paths to the ray files given a step
    # rays_list = os.listdir(rays_path)[::step]
    rays_list = np.array(os.listdir(rays_path))[[0, 2, 4, 7, 11, 15, 19, 23, 28]]
    # How many files are we considering given a step?
    rays_list_len = len(rays_list)

    if key == 0:
        flag = 4  # The desired number of subplots on one Figure
    elif key == 1:
        flag = 9
    else:
        sys.exit("key argument must be equal to 0 or 1!")

    # What is the remainder of dividing the number of files with rays by flag?
    dr = rays_list_len % flag

    # List of column names that are contained within each raysXX file
    columns = ['lon1', 'lat1', 'z1', 'lon2', 'lat2', 'z2', 'Vg']

    # Scale limits for colorbar
    # A class which, when called, linearly normalizes data into the
    # ``[0.0, 1.0]`` interval.
    norm = plt.Normalize(vmin=-60, vmax=60)

    for i in range(0, rays_list_len):
        # Load data and use them to create a DataFrame object
        cur_filename = rays_list[i]
        cur_period = int(cur_filename[4:6])
        cur_data = np.loadtxt(fname=rays_path+'\\'+cur_filename)
        cur_rays = pd.DataFrame(data=cur_data, columns=columns)

        # Calculate the percentage of deviation of each group velocity
        # from the average group velocity for the current period
        cur_rays['Deviaton from aVg, %'] = 100 * \
            (cur_rays['Vg'] - cur_rays['Vg'].mean()) / cur_rays['Vg'].mean()

        # Create a suitable Figure object
        # When it is impossible to place 4 subplots on a Figure
        if (dr != 0) & (i == (rays_list_len - dr)):
            if dr == 1:
                sys.exit("Drawing a Figure with a single subplot isn't supported")

            fig, ax = plt.subplots(1, dr, dpi=150, figsize=(16, 10))

            if dont_open:
                plt.close(fig)  # So that the figure is not shown

            # Make 1d vector from 2d ax array
            ax = ax.flatten()
            # Axis labels
            ax[0].set(xlabel='longitude, degrees', ylabel='latitude, degrees')
            for q in range(1, dr):
                ax[q].set(xlabel='longitude, degrees', yticklabels=[])
            # Adjust the distance between subplots
            fig.subplots_adjust(wspace=0.08)
            # Colorbar setup
            cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=plt.cm.RdBu),
                              label='deviation from the average group velocity, %',
                              orientation='horizontal', aspect=50, ax=ax, pad=0.1)

            # Reset the flag
            flag = 0

        # When it is possible to place 4 subplots on a Figure
        elif (key == 0) & (flag == 4):
            fig, ax = plt.subplots(2, 2, dpi=150, figsize=(16, 10))

            if dont_open:
                plt.close(fig)  # So that the figure is not shown

            # Make 1d vector from 2d ax array
            ax = ax.flatten()
            # Axis labels
            ax[0].set(ylabel='latitude, degrees', xticklabels=[])
            ax[1].set(xticklabels=[], yticklabels=[])
            ax[2].set(xlabel='longitude, degrees', ylabel='latitude, degrees')
            ax[3].set(xlabel='longitude, degrees', yticklabels=[])
            # Adjust the distance between subplots
            fig.subplots_adjust(wspace=0.08, hspace=0.045)
            # Colorbar setup
            fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=plt.cm.RdBu),
                         label='deviation from the average group velocity, %',
                         aspect=30, ax=ax)

            # Reset the flag
            flag = 0

        elif (key == 1) & (flag == 9):
            fig, ax = plt.subplots(3, 3, dpi=150, figsize=(16, 10))

            if dont_open:
                plt.close(fig)  # So that the figure is not shown

            # Make 1d vector from 2d ax array
            ax = ax.flatten()
            # Axis labels
            ax[0].set(ylabel='latitude, degrees', xticklabels=[])
            ax[1].set(xticklabels=[], yticklabels=[])
            ax[2].set(xticklabels=[], yticklabels=[])
            ax[3].set(ylabel='latitude, degrees', xticklabels=[])
            ax[4].set(xticklabels=[], yticklabels=[])
            ax[5].set(xticklabels=[], yticklabels=[])
            ax[6].set(xlabel='longitude, degrees', ylabel='latitude, degrees')
            ax[7].set(xlabel='longitude, degrees', yticklabels=[])
            ax[8].set(xlabel='longitude, degrees', yticklabels=[])
            # Adjust the distance between subplots
            fig.subplots_adjust(wspace=0.08, hspace=0.045)
            # Colorbar setup
            fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=plt.cm.RdBu),
                         label='deviation from the average group velocity, %',
                         aspect=30, ax=ax)

            # Reset the flag
            flag = 0

        # After creating a suitable Figure, we will loop through
        # each subplot on it using flag

        # Set aspect ratio of the current axis
        ax[flag].set_aspect(aspect=yx_aspect_ratio)
        # Set axis limits
        ax[flag].set(xlim=x_limits, ylim=y_limits)
        # Draw ticks on all 4 axes
        ax[flag].tick_params(labelsize=8, top=True, right=True)

        # Plot rays
        for j in range(cur_rays.shape[0]):
            x = [cur_rays.loc[j, 'lon1'], cur_rays.loc[j, 'lon2']]
            y = [cur_rays.loc[j, 'lat1'], cur_rays.loc[j, 'lat2']]
            ax[flag].plot(x, y, linewidth=0.75,
                          color=plt.cm.RdBu(norm(cur_rays.loc[j, 'Deviaton from aVg, %'])),
                          alpha=1, zorder=j)
        # Plot coastal line
        for m in range(len(coastlines_list)):
            cur_line = coastlines_list[m]
            ax[flag].plot(cur_line['longitude'], cur_line['latitude'],
                          color='black', linewidth=1.7, zorder=-1)
        # Plot stations
        for k in range(stations.shape[0]):
            ax[flag].scatter(stations.iloc[k, 0], stations.iloc[k, 1],
                             marker='^', s=50, color='white', edgecolor='black',
                             zorder=j+1)
        # Annotate the current plot
        if key == 0:
            ax[flag].text(0.87, 0.05, f"period: {cur_period} s",
                          horizontalalignment='center', verticalalignment='center',
                          transform=ax[flag].transAxes)
        elif key == 1:
            ax[flag].text(0.84, 0.05, f"period: {cur_period} s", size='small',
                          horizontalalignment='center', verticalalignment='center',
                          transform=ax[flag].transAxes)

        flag += 1

        # Save Figure if required
        if ((save_path != 0)
            & (((flag == 4) & (key == 0))
                | ((flag == 9) & (key == 1))
                | ((dr != 0) & (i == (rays_list_len - 1))))):
            iterator = 1
            while True:
                if os.path.exists(save_path + f"\\rays_{iterator}.png"):
                    iterator += 1
                else:
                    break
            output_file = save_path + f"\\rays_{iterator}.png"
            fig.savefig(output_file, bbox_inches='tight', dpi=150)
            print("The current figure saved to " + output_file)


def visualize_velocity_deviations_ru(coastlines_list, x_limits, y_limits,
                                  stations_path, rays_path,
                                  step=1, yx_aspect_ratio=1, key=0,
                                  dont_open=0, save_path=0, ru=False):
    '''
    For a given set of periods, creates maps with coastlines, seismic stations
    and rays between them.
    The rays are colored according to the Rayleigh group velocity.
    NOTE: the behavior when one subplot is placed on a Figure is unscripted;
    all captions will be made in Russian.

    Parameters:

    coastlines_list - A list of coastlines
    that lie within a previously specified range of geographic coordinates

    x_limits, y_limits - the geographic coordinates within which maps are plotted

    stations_path - the path to the excel-file with information
    about the relevant stations

    rays_path - the path to the directory containing the ray files

    step - visualize every step-file in rays_path

    yx_aspect_ratio - the ratio of the length of the unit segment of the y-axis
    to the unit segment of the x-axis;
    by default, these segments are equal in length

    key - 0 for Figures with 2x2 subplots, 1 for Figures with 3x3 subplots

    dont_open - do not open or open windows with built visualizations;
    by default, the windows will open

    save_path - the path to the directory where Figures will be saved;
    by default, Figures are not saved
    '''

    # Create DataFrame with information about the relevant stations
    stations = pd.read_excel(stations_path, index_col=0)

    # Check if step is legit
    if type(step) != int:
        sys.exit("step argument must be an integer!")
    # Create a list with the paths to the ray files given a step
    rays_list = os.listdir(rays_path)[::step]
    # How many files are we considering given a step?
    rays_list_len = len(rays_list)

    if key == 0:
        flag = 4  # The desired number of subplots on one Figure
    elif key == 1:
        flag = 9
    else:
        sys.exit("key argument must be equal to 0 or 1!")

    # What is the remainder of dividing the number of files with rays by flag?
    dr = rays_list_len % flag

    # List of column names that are contained within each raysXX file
    columns = ['lon1', 'lat1', 'z1', 'lon2', 'lat2', 'z2', 'Vg']

    # Scale limits for colorbar
    # A class which, when called, linearly normalizes data into the
    # ``[0.0, 1.0]`` interval.
    norm = plt.Normalize(vmin=-60, vmax=60)

    for i in range(0, rays_list_len):
        # Load data and use them to create a DataFrame object
        cur_filename = rays_list[i]
        cur_period = int(cur_filename[4:6])
        cur_data = np.loadtxt(fname=rays_path+'\\'+cur_filename)
        cur_rays = pd.DataFrame(data=cur_data, columns=columns)

        # Calculate the percentage of deviation of each group velocity
        # from the average group velocity for the current period
        cur_rays['Deviaton from aVg, %'] = 100 * \
            (cur_rays['Vg'] - cur_rays['Vg'].mean()) / cur_rays['Vg'].mean()

        # Create a suitable Figure object
        # When it is impossible to place 4 subplots on a Figure
        if (dr != 0) & (i == (rays_list_len - dr)):
            if dr == 1:
                sys.exit("Drawing a Figure with a single subplot isn't supported")

            fig, ax = plt.subplots(1, dr, dpi=150, figsize=(16, 10))

            if dont_open:
                plt.close(fig)  # So that the figure is not shown

            # Make 1d vector from 2d ax array
            ax = ax.flatten()
            # Axis labels
            ax[0].set(xlabel='долгота, градусы', ylabel='широта, градусы')
            for q in range(1, dr):
                ax[q].set(xlabel='долгота, градусы', yticklabels=[])
            # Adjust the distance between subplots
            fig.subplots_adjust(wspace=0.08)
            # Colorbar setup
            cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=plt.cm.RdBu),
                              label='отклонение от средней групповой скорости, %',
                              orientation='horizontal', aspect=50, ax=ax, pad=0.1)

            # Reset the flag
            flag = 0

        # When it is possible to place 4 subplots on a Figure
        elif (key == 0) & (flag == 4):
            fig, ax = plt.subplots(2, 2, dpi=150, figsize=(16, 10))

            if dont_open:
                plt.close(fig)  # So that the figure is not shown

            # Make 1d vector from 2d ax array
            ax = ax.flatten()
            # Axis labels
            ax[0].set(ylabel='широта, градусы', xticklabels=[])
            ax[1].set(xticklabels=[], yticklabels=[])
            ax[2].set(xlabel='долгота, градусы', ylabel='широта, градусы')
            ax[3].set(xlabel='долгота, градусы', yticklabels=[])
            # Adjust the distance between subplots
            fig.subplots_adjust(wspace=0.08, hspace=0.045)
            # Colorbar setup
            fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=plt.cm.RdBu),
                         label='отклонение от средней групповой скорости, %',
                         aspect=30, ax=ax)

            # Reset the flag
            flag = 0

        elif (key == 1) & (flag == 9):
            fig, ax = plt.subplots(3, 3, dpi=150, figsize=(16, 10))

            if dont_open:
                plt.close(fig)  # So that the figure is not shown

            # Make 1d vector from 2d ax array
            ax = ax.flatten()
            # Axis labels
            ax[0].set(ylabel='широта, градусы', xticklabels=[])
            ax[1].set(xticklabels=[], yticklabels=[])
            ax[2].set(xticklabels=[], yticklabels=[])
            ax[3].set(ylabel='широта, градусы', xticklabels=[])
            ax[4].set(xticklabels=[], yticklabels=[])
            ax[5].set(xticklabels=[], yticklabels=[])
            ax[6].set(xlabel='долгота, градусы', ylabel='широта, градусы')
            ax[7].set(xlabel='долгота, градусы', yticklabels=[])
            ax[8].set(xlabel='долгота, градусы', yticklabels=[])
            # Adjust the distance between subplots
            fig.subplots_adjust(wspace=0.08, hspace=0.045)
            # Colorbar setup
            fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=plt.cm.RdBu),
                         label='отклонение от средней групповой скорости, %',
                         aspect=30, ax=ax)

            # Reset the flag
            flag = 0

        # After creating a suitable Figure, we will loop through
        # each subplot on it using flag

        # Set aspect ratio of the current axis
        ax[flag].set_aspect(aspect=yx_aspect_ratio)
        # Set axis limits
        ax[flag].set(xlim=x_limits, ylim=y_limits)
        # Draw ticks on all 4 axes
        ax[flag].tick_params(labelsize=8, top=True, right=True)

        # Plot rays
        for j in range(cur_rays.shape[0]):
            x = [cur_rays.loc[j, 'lon1'], cur_rays.loc[j, 'lon2']]
            y = [cur_rays.loc[j, 'lat1'], cur_rays.loc[j, 'lat2']]
            ax[flag].plot(x, y, linewidth=0.75,
                          color=plt.cm.RdBu(norm(cur_rays.loc[j, 'Deviaton from aVg, %'])),
                          alpha=1, zorder=j)
        # Plot coastal line
        for m in range(len(coastlines_list)):
            cur_line = coastlines_list[m]
            ax[flag].plot(cur_line['longitude'], cur_line['latitude'],
                          color='black', linewidth=1.7, zorder=j+1)
        # Plot stations
        for k in range(stations.shape[0]):
            ax[flag].scatter(stations.iloc[k, 0], stations.iloc[k, 1],
                             marker='^', s=50, color='white', edgecolor='black',
                             zorder=j+2)
        # Annotate the current plot
        if key == 0:
            ax[flag].text(0.87, 0.05, f"период: {cur_period} с",
                          horizontalalignment='center',
                          verticalalignment='center',
                          transform=ax[flag].transAxes)
        elif key == 1:
            ax[flag].text(0.84, 0.05, f"период: {cur_period} с",
                          size='small',
                          horizontalalignment='center',
                          verticalalignment='center',
                          transform=ax[flag].transAxes)

        flag += 1

        # Save Figure if required
        if ((save_path != 0)
            & (((flag == 4) & (key == 0))
                | ((flag == 9) & (key == 1))
                | ((dr != 0) & (i == (rays_list_len - 1))))):
            iterator = 1
            while True:
                if os.path.exists(save_path + f"\\rays_{iterator}.png"):
                    iterator += 1
                else:
                    break
            output_file = save_path + f"\\rays_{iterator}.png"
            fig.savefig(output_file, bbox_inches='tight', dpi=150)
            print("The current figure saved to " + output_file)

def extract_coastlines(coastlines_path, x_min, x_max, y_min, y_max):
    '''
    Сreates a list of coastlines
    that lie in a given range of geographic coordinates.

    Parameters:

    path_coastlines - path to a .bln file
    that contains information about all the coastlines of the planet

    x_min, x_max, y_min, y_max - the geographic coordinates
    within which coastline information is extracted and stored
    '''

    # Load coastlines for the entire planet
    coastlines = pd.read_csv(filepath_or_buffer=coastlines_path, header=None,
                             names=['longitude', 'latitude'], sep="\s+")

    # Each coastline in coastlines DataFrame starts with a line that contains
    # the number of points and the missing value.
    # We want to create an object whose elements contain
    # indexes of the coastlines file corresponding to missing values.
    nan_locations = coastlines.isna().sum(axis=1)  # pd.Series with values 0 or 1
    nan_locations = nan_locations[nan_locations == 1].index.values # np.array

    coastlines_list = []  # list for storing DataFrames with individual lines

    # We will look at the chunks of the coastlines file between the NaN-values.
    # These chunks correspond to the individual lines.
    beginning = 0   # first index of the current chunk
    end = nan_locations[1]  # last index of the current chunk
    flag = 1  # iterator for nan_locations

    while True:
        # the last value of a slice is end - 1!
        cur_line = coastlines.iloc[beginning:end]
        # We consider only those values that lie within the given limits
        cur_line = cur_line[(cur_line['longitude'] >= x_min)
                            & (cur_line['longitude'] <= x_max)
                            & (cur_line['latitude'] >= y_min)
                            & (cur_line['latitude'] <= y_max)]

        # Check if there are any valid values left in cur_line
        if not cur_line.empty:
            print("\nSomething was found!\nAppending to the list...")
            # If something are still in cur_line, append it to coastlines_list
            coastlines_list.append(cur_line)

        # Iterate until flag is equal to (the last index of nan_locations) + 1
        # Remember: the first index is 0, not 1!
        flag += 1
        if flag == len(nan_locations):
            break

        # Redefine the beginning and end indices of the coastlines slice
        beginning = end + 1
        end = nan_locations[flag]

    return coastlines_list


def visualize_wave_propagation(crosscorr_matrix_path, distances_path,
                               srate, frange, lpbutter_order, hpbutter_order,
                               amp_coeff=1, alpha=0.7, fill_color='black',
                               xlim=None, ylim=None, dont_open=0, save_path=0,
                               ru=False):
    '''
    Creates a figure with filtered cross-correlation functions
    plotted depending on the distance between pairs of stations.

    Parameters:

    crosscorr_matrix_path - full path to the DataFrame with cross-correlations

    distances_path - full path to the DataFrame with distances
    between all station pairs

    srate - sampling rate of the cross-correlation data

    frange - frequency band of the bandpass (lowpass + highpass Butter) filter
    that will be applied to the data before their visualization

    lpbutter_order, hpbutter_order - orders of Butterworth low- and high-pass filters,
    which will be applied to the data before their visualization

    amp_coeff - cross-correlation amplitude scaling factor;
    the default value is 1

    alpha - the value used for blending;
    alpha must be in the range 0-1, inclusive.

    fill_color - what color to use to fill areas
    where the cross-correlation amplitude is positive;
    the default value is 'black'

    xlim, ylim -sets the limits of the x and y axes; if none are given,
    the function uses xlim=[-150, 150], ylim=[-10, 1.05*max_dist].

    dont_open - open (0) or don't open (1) the window with created visualization;
    by default, the window will open

    save_path - the path to the directory in which Figure will be saved;
    by default, Figure is not saved

    ru - if true, all captions will be made in Russian
    '''

    # Create DataFrames with cross-correlations
    # and distances between all station pairs
    crosscorr_matrix = pd.read_excel(io=crosscorr_matrix_path, index_col=0)
    t = crosscorr_matrix.index.values  # time vector
    distances = pd.read_excel(io=distances_path, index_col=0)

    # Create filter coefficients. They can be applied to all columns from our data!
    lp_filter = signal.butter(N=lpbutter_order, Wn=frange[1], btype='lowpass',
                              fs=srate, output='sos')
    hp_filter = signal.butter(N=hpbutter_order, Wn=frange[0], btype='highpass',
                              fs=srate, output='sos')

    # The maximum distance between stations in station pairs
    max_dist = distances.values.max()/1000
    # print(max_dist)

    fig, ax = plt.subplots(dpi=200, figsize=(16, 9))

    if dont_open:
        plt.close(fig)  # So that the figure is not shown

    for column in crosscorr_matrix.columns:
        # Filter current column
        y = signal.sosfiltfilt(lp_filter, crosscorr_matrix[column].values)
        y = signal.sosfiltfilt(hp_filter, y)
        # Normalization of the result
        # Z-score multiplied by additional scaling coeff
        y = amp_coeff*stats.zscore(y)
        # Distance between the current pair of stations
        dist = distances.loc[column].values/1000
        y += dist
        # Plot a line and then fill in its parts that are larger than 0
        ax.plot(t, y, color='black', linewidth=0.5, alpha=alpha)
        ax.fill_between(x=t, y1=y, y2=dist, where=y > dist,
                        color=fill_color, edgecolor=None,
                        alpha=alpha)

    # If xlim and ylim are not defined, use these default axis limits
    if (not xlim) & (not ylim):
        if ru:
            ax.set(xlim=[-150, 150], ylim=[-10, 1.05*max_dist],
                   xlabel='время, с',
                   ylabel='расстояние между станциями в паре, км')
        else:
            ax.set(xlim=[-150, 150], ylim=[-10, 1.05*max_dist],
                   xlabel='time, s',
                   ylabel='distance between station pairs, km')
    # If xlim and ylim are incorrectly defined by the user, show this error msg
    elif (len(xlim) != 2) or (len(ylim) != 2):
        sys.exit("Something is wrong with xlim / ylim")
    # Otherwise, use these user defined axis limits
    else:
        if ru:
            ax.set(xlim=xlim, ylim=ylim,
                   xlabel='время, с',
                   ylabel='расстояние между станциями в паре, км')
        else:
            ax.set(xlim=xlim, ylim=ylim,
                   xlabel='time, s',
                   ylabel='distance between station pairs, km')

    # Save Figure if required
    if save_path != 0:
        output_file = save_path + "surface_wave_propagation.png"
        fig.savefig(output_file, bbox_inches='tight', dpi=300)
        print("Figure is saved to " + output_file)
