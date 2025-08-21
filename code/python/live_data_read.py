import serial
from time import time
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
import numpy as np

ser = serial.Serial('COM7', baudrate=115200)
print(ser.name)

x = []
t = []
I = []

num_peaks = 2
times = 0

ser.readline() # read and do nothing with first line
start = time()
while True:
    line = ser.readline().decode('utf-8').strip()  # Read a line and decode it
    print(line)
    if not line:  # Stop if no data is received
        break

    x.append(float(line.split()[0])/1000) # convert mm to m
    t.append((time() - start))
    I.append(float(line.split()[2])/1000) # convert from mA to A
    # print(x)
    # print(len(x))
    # print(len(x)%10)
    if len(x) >= 100 and len(x)%10 == 0:
        x100 = x[-100:]
        t100 = t[-100:]
        current = np.mean(I[-100:])
        
        # print(np.round(x100, 3), len(x100))
        # print(np.round(t100, 3), len(t100))
        # print(current)

        amplitudes, frequencies = extracter(x100, t100, num_peaks)
        # ⇓ replace with constitutive relationship between damping and current ⇓
        damping = 1.5*(current**2)
        # ⇑--------------------------------------------------------------------⇑
        power = 0
        for i in range(len(amplitudes)):
            power += (1/2)*damping*(amplitudes[i]**2)*((frequencies[i]*2*np.pi)**2)
        print(f"\nPower: {round(power, 4)}W")
        print(f"Peaks: amplitude of {np.round(amplitudes, 4)}m and frequency of {np.round(frequencies, 4)}Hz; Mean current: {current}A\n")