# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 17:45:46 2022

@author: Igor
"""
import numpy as np
import pandas as pd
from scipy import signal
from functions_visualization import visualize_iir_filter, visualize_wave_propagation

#
crosscorr_matrix_path = "E:/Work/Projects/Central_Kamchatka/Cross-correlations/averaged_crosscorrs_full.xlsx"
distances_path = "E:/Work/Projects/Central_Kamchatka/Temp/stations_distance.xlsx"
save_path = 'E://Work//Projects//Central_Kamchatka//_Paper//Figures//'

srate = 10  # the data's sampling rate
frange = [0.05, 0.1]  # desired frequency range
# %% Search for a Butterworth filter that is applicable
# to the available data

# Define some variables for calculations and plotting of filters
crosscorr_matrix = pd.read_excel(io=crosscorr_matrix_path, index_col=0)
npts = crosscorr_matrix.shape[0]
srate = 10
hz = np.linspace(0, srate/2, int(np.floor(npts/2 + 1)))
# Desired frequency range
frange = [0.05, 0.1]

# Ideal shape of the highpass filter
hpfilter_xshape = [0, frange[0], frange[0], srate/2]
hpfilter_yshape = [0, 0, 1, 1]
# Ideal shape of the lowpass filter
lpfilter_xshape = [0, frange[1], frange[1], srate/2]
lpfilter_yshape = [1, 1, 0, 0]
freq_xlim = [0, 1.1]

# Highpass Butter's order selection
for i in range(2, 11):
    # Current order of the Butterworth filter
    hpbutter_order = i
    # Create filter coefficients
    hpbutter_b, hpbutter_a = signal.butter(hpbutter_order, frange[0],
                                           btype='highpass', fs=srate)
    # Plot
    visualize_iir_filter(hpbutter_b, hpbutter_a, 'highpass Butterworth filter',
                         hpbutter_order, 0, hpfilter_xshape, hpfilter_yshape,
                         npts, hz, freq_xlim=[0, 0.5])

# Lowpass Butter's order selection
for i in range(5, 16):
    # Current order of the Butterworth filter
    lpbutter_order = i
    # Create filter coefficients
    lpbutter_b, lpbutter_a = signal.butter(lpbutter_order, frange[1],
                                           btype='lowpass', fs=srate)
    # Plot
    visualize_iir_filter(lpbutter_b, lpbutter_a, 'lowpass Butterworth filter',
                         lpbutter_order, 1, lpfilter_xshape, lpfilter_yshape,
                         npts, hz, freq_xlim=[0, 0.5])
# %%
# Chosen Butterworth's filter orders
lpbutter_order = 8
hpbutter_order = 7
# %%
visualize_wave_propagation(crosscorr_matrix_path, distances_path, srate, frange,
                           lpbutter_order, hpbutter_order,
                           amp_coeff=0.75, fill_color='black',
                           save_path=save_path, ru=True)
