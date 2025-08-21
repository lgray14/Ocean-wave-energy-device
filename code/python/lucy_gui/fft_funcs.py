import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
from numpy.fft import fft, fftfreq
from scipy.signal import windows, find_peaks

def quadratic_peak_interpolation(amps, freqs, idx):
    """
    Refine peak frequency and amplitude using quadratic interpolation.

    Parameters:
        amps (ndarray): FFT amplitude spectrum
        freqs (ndarray): Corresponding frequency bins
        idx (int): Index of the peak

    Returns:
        freq_peak (float): Interpolated frequency
        amp_peak (float): Interpolated amplitude
        p (float): Offset in bins from center
    """
    p = 0.5 * (amps[idx - 1] - amps[idx + 1]) / (amps[idx - 1] - 2*amps[idx] + amps[idx + 1])
    freq_peak = freqs[idx] + p*(freqs[1] - freqs[0])
    amp_peak = amps[idx] - 0.25*(amps[idx - 1] - amps[idx + 1])*p
    return freq_peak, amp_peak

def extracter(x, t, k):
    """
    Extract k peak ampliudes and their corresponding frequencies using zero-padded fast Fourier transform.

    Parameters:
        x: original displacement signal
        t: corresponding time data
        k: number of peak frequencies to be identified

    Returns:
        peak_amps: identified peak amplitudes, list of length k
        peak_freqs: corresponding frequencies
    """

    n = len(t)
    window = windows.hann(n)
    y_windowed = x*window

    # Zero-padding
    zero_pad_factor = 8
    n_padded = n*zero_pad_factor

    # do fft on zero-padded data
    fft_vals = fft(y_windowed, n=n_padded)
    freqs_raw = fftfreq(n_padded, d=t[1] - t[0])
    # take the second half/valid data
    freqs = freqs_raw[:n_padded//2]
    amps = np.abs(fft_vals[:n_padded//2]) / (np.sum(window)/2)

    # find peaks to analyze around
    peaks, _ = find_peaks(amps)
    peak_freqs = []
    peak_amps = []

    combined_lists = list(zip(amps[peaks], freqs[peaks], peaks))
    sorted_lists = sorted(combined_lists, key=lambda x: x[0]) # sort frequencies and amplitudes AND peak indices based on amplitude

    # 1. list sorts by amplitude from smallest to largest
    # 2. the peak will always be the third item in the row
    # therefore, the index of the largest peak should be sorted_lists[-1][2]
    # ==> to index in for loop, next to add 1 to i to correct for zero and should be the same -- sorted_lists[-(i+1)][2]
    # i rows from the bottom, third entry in row

    for i in range(k):
    # Perform quadratic interpolation
        freq_peak, amp_peak = quadratic_peak_interpolation(amps, freqs, sorted_lists[-(i+1)][2]) # get the peak index for the 
        peak_freqs.append(freq_peak)
        peak_amps.append(amp_peak)
    
    return peak_amps, peak_freqs