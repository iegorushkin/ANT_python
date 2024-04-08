# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 15:49:29 2021

@author: Igor
"""
import numpy as np
import copy
from scipy import signal, fftpack
from obspy.signal.invsim import cosine_taper
from copy import deepcopy


def passband_butter(st, frange, hp_order, lp_order):
    '''
    Passband filters a single trace from the obspy stream object st
    using high- and lowpass Butterworth filters.

    Parameters:

    st - obspy stream object containing only 1 trace

    frange - сombined filter bandwidth

    hp_order - highpass Butterworth filter order

    lp_order - lowpass Butterworth filter order

    Returns an obspy stream object with ONE filtered trace.
    '''

    # Define variables for convenience
    sig = st[0].data
    srate = int(st[0].stats.sampling_rate)

    # Lets filter!
    # create coefficients for 2 Butterworth filters - high- and lowpass
    hp_sos = signal.butter(hp_order, frange[0], btype='highpass',
                           output='sos', fs=srate)
    lp_sos = signal.butter(lp_order, frange[1], btype='lowpass',
                           output='sos', fs=srate)

    # Filtering the signal with one filter, then the other.
    sig = signal.sosfiltfilt(hp_sos, sig)
    sig = signal.sosfiltfilt(lp_sos, sig)

    # Save the filtered signal in the obspy stream object st
    st[0].data = sig

    return st


def spect_norm(st):
    '''
    Computes the spectral normalization of the first trace
    of an obspy stream object.
    This is done by dividing the data spectrum by its absolute value.

    Parameters:

    st - obspy stream object containing only 1 trace

    Returns an obspy stream object with ONE spectrally normalized trace.
    '''

    # Define variables for convenience
    sig = st[0].data

    # taper data
    sig *= cosine_taper(len(sig), 0.01)

    # Fast Fourier Transform
    sig_fft = fftpack.fft(sig)
    sig_fft_norm = sig_fft / np.abs(sig_fft)

    # Inverse Fast Fourier Transform
    sig = np.real(fftpack.ifft(sig_fft_norm))
    st[0].data = sig

    return st


def spect_whiten(st, frange, file_name='', npts_multiplier=False):
    '''
    Flattens the seismogram spectrum in the specified frequency range
    and zeroes out outside it.
    Introduces a short smooth transition zone between
    the normalized and zeroed frequency ranges to avoid time domain artifacts.

    Parameters:

    st - obspy stream object containing only 1 trace

    frange - spectrum normalization boundaries

    file_name - name of the mseed file from which the stream is read

    npts_multiplier - whether the whitened spectrum should be multiplied
    by the number of points in the signal.

    Returns an obspy stream object with ONE trace
    that is spectrally normalized within the specified boundaries.
    '''

    # Define variables for convenience
    sig = st[0].data
    srate = int(st[0].stats.sampling_rate)
    npts = len(sig)
    hz = np.linspace(0, srate/2, int(np.floor(npts/2) + 1))

    # Frequency range where whitening needs to be applied
    freqs = frange
    # Indices corresponding to this frequency range.
    ifreqs = np.where((hz >= freqs[0]) & (hz <= freqs[1]))[0]
    # Indices of min and max frequencies of the desired range
    imin = ifreqs[0]
    imax = ifreqs[-1]

    # Taper lenght
    npts_taper = 120

    # Start of the left transition zone
    l_trans = imin - npts_taper
    # Check boundries
    if l_trans <= 0:
        l_trans = 1

    # End of the right transition zone
    r_trans = imax + npts_taper
    # Check boundries
    if r_trans >= int(np.floor(npts/2) + 1):
        r_trans = int(np.floor(npts/2) + 1)

    # Signal FFT
    fft_sig = fftpack.fft(sig)

    # Zeros before l_trans and after r_trans
    fft_sig[:l_trans] *= 0
    fft_sig[r_trans+1:] *= 0

    if npts_multiplier:
        # Left tapering
        fft_sig[l_trans:imin] = (
            npts*np.sin(np.linspace(0, np.pi/2, imin - l_trans))**2
            * np.exp(1j*np.angle(fft_sig[l_trans:imin]))
              )
        # Right tapering
        fft_sig[imax+1:(r_trans+1)] = (
            npts*np.cos(np.linspace(0, np.pi/2, r_trans - imax))**2
            * np.exp(1j*np.angle(fft_sig[imax+1:(r_trans+1)]))
            )
        # Pass band
        fft_sig[imin:imax+1] = npts*np.exp(1j*np.angle(fft_sig[imin:imax+1]))
    else:
        # Left tapering
        fft_sig[l_trans:imin] = (
            np.sin(np.linspace(0, np.pi/2, imin - l_trans))**2
            * np.exp(1j*np.angle(fft_sig[l_trans:imin]))
              )
        # Right tapering
        fft_sig[imax+1:(r_trans+1)] = (
            np.cos(np.linspace(0, np.pi/2, r_trans - imax))**2
            * np.exp(1j*np.angle(fft_sig[imax+1:(r_trans+1)]))
            )
        # Pass band
        fft_sig[imin:imax+1] = np.exp(1j*np.angle(fft_sig[imin:imax+1]))

    # Hermitian symmetry
    fft_sig[-len(hz)+1::] = np.conjugate(fft_sig[1:len(hz)])[::-1]
    # Check if the values really are symmetrical around Nyquist.
    try:
        print("\nChecking Hermitian symmetry after whitening")
        print(f"File: {file_name}")
        print(f'Nyquist + 16333: {fft_sig[len(hz)-1+16333]}')
        print(f'Nyquist - 16333: {fft_sig[len(hz)-1-16333]}')
    except IndexError:
        print("\nTrace is too short to test the Hermitian symmetry.")
        print(f"{len(fft_sig)}")

    # For testing of taper
    print(f"\nLeft taper length: {imin - l_trans},"
          + f" right taper length: {r_trans - imax}")

    # Inverse FFT
    sig = np.real(fftpack.ifft(fft_sig))
    st[0].data = sig

    return st


def onebit_norm(st):
    '''
    Replaces the data stored in a single trace of an obspy stream object st
    with its element-wise sign indication.

    Parameters:

    st - obspy stream object containing only 1 trace

    Returns an obspy stream object with its single data trace replaced with
    -1 if x < 0, 0 if x==0, 1 if x > 0.
    '''
    # Define variables for convenience
    sig = st[0].data

    # one-bit normalization
    sig = np.sign(sig)
    st[0].data = sig

    return st

def ram_norm(st, max_period):
    """
    Performs time-domain data normalization based on the running absolute mean
    weighting scheme.
    Important note: it is assumed that each weight is calculated independently.

    Parameters:

    st - obspy stream object containing only 1 trace

    max_period - the minimum frequency present in the (filtered) data;
    used to calculate the length of a running window.

    Returns an obspy stream object with its single data trace normalized.
    """
    # Define variable(s) for convenience
    srate = st[0].stats.sampling_rate  # sampling rate
    max_period_length = int(max_period*srate)  # 1/2 of windows size in samples
    # Variables for the signal and its length
    norm_sig = deepcopy(st[0].data)
    norm_sig_length = len(norm_sig)
    # Array for storing the calculated weights.
    weights = np.zeros(norm_sig_length)

    # Calculation of weights for all signal values and storing them to an array
    for i in range(norm_sig_length):
        # The calculation of a moving window is performed
        # taking into account that the point with i-index
        # should be located in its center.
        # However, at the edges of the data, in order to avoid the loss
        # of useful information, special conditions are set.
        lower_bound = max(0, i - max_period_length)
        upper_bound = min(norm_sig_length - 1, i + max_period_length)
        weights[i] = np.mean(np.abs(norm_sig[lower_bound:(upper_bound + 1)]))

    # Dividing the signal points by the corresponding weights.
    norm_sig /= weights

    st[0].data = norm_sig
    return st
