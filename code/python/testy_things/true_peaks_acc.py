import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
from numpy.random import rand
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

sr = 50 #Hz -- sampling rate
# sampling interval
ts = 1.0/sr # timestamp
time_range = 10
t = np.arange(0,time_range,ts)

accuracies = []
for i in range(100):
    # Generate fake signal
    signals = 1
    x = 0
    real_power = 0
    for i in range(signals):
        frequency = round(sr/8*rand(), 3)
        amp = round(rand(), 2)
        x += amp*np.sin(2*np.pi*frequency*t)
        # print("Freq %s:" %i, frequency, "Amp %s:" %i, amp)
        real_power += 1/2*((frequency*2*np.pi)**2)*(amp**2)


    # READ GUI RESULTS TRIAL
    # x = []
    # t = []
    # with open("gui_results/gui_real_test2.txt", "r") as f:
    #     lines = f.readlines()

    #     # Loop through all lines, ignoring header.
    #     # Add last element to list (i.e. the process name)
    #     for l in lines[2:]:
    #         x.append((float(l.split()[1]))/100) # take cm data and covert to m
    #         t.append(float(l.split()[0]))

    # # find average time gap ==> sampling frequency
    # diffs = []
    # for i in range(len(t)-1):
    #     diffs.append(t[i+1] - t[i])
    # avg_gap = np.mean(diffs)
    # print(avg_gap)
    # END

    # plt.figure(figsize = (9, 5))
    # plt.subplot(311)
    # plt.plot(t, x, label='Original Signal')
    # plt.title('Original Signal')
    # plt.ylabel('Displacement (m)')

    n = len(t)

    # Window the signal
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

    # ---- PLOT FFT RESULTS ----
    # plt.subplot(312)
    # plt.stem(freqs, amps, markerfmt=" ", basefmt="-b", label='Zero-padded FFT')
    # # plt.axvline(freq_peak, color='r', linestyle='--', label=f'Interp. Peak ≈ {freq_peak:.4f} Hz')
    # plt.title("Frequency Spectrum from FFT")
    # plt.xlabel("Frequency (Hz)")
    # plt.ylabel("Amplitude")
    # plt.legend()

    #---------------- just one highest peak ----------------
    # # Find the index of the max peak
    # idx = np.argmax(amps)
    # if idx <= 0 or idx >= len(amps) - 1:
    #     raise ValueError("Peak is at FFT boundary; cannot interpolate.")
    # freq_peak, amp_peak = quadratic_peak_interpolation(amps, freqs, idx)

    # reconstructed = amp_peak*np.sin(2*np.pi*freq_peak*t) # reconstructed signal
    # # --- PRINT RESULTS ---
    # print(f"Detected frequency: {freq_peak:.4f} Hz")
    # print(f"Estimated amplitude: {amp_peak:.4f}")
    #-------------------------------------------------------

    # ~~~~~~~~~~~~ k highest peaks ~~~~~~~~~~~~
    k = 1 # filter to k highest peaks
    peaks, _ = find_peaks(amps)
    # print(amps)
    # print(peaks)
    # print(np.round(amps[peaks], 2))
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
    # reconstructed singal
    reconstructed = 0
    modelled_power = 0
    for i in range(len(peak_freqs)):
        reconstructed += peak_amps[i]*np.sin(2*np.pi*peak_freqs[i]*t)
        modelled_power += 1/2*((peak_freqs[i]*2*np.pi)**2)*(peak_amps[i]**2)

    # --- PRINT RESULTS ---
    # print(f"Detected frequency: {np.round(peak_freqs, 4)} Hz")
    # print(f"Estimated amplitude: {np.round(peak_amps, 4)} m")
    # print(f"Real power: {round(real_power, 4)}W/damping vs Modelled power: {round(modelled_power, 4)}W/damping")
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    accuracies.append(modelled_power/real_power)
    # # ---- PLOT RECONSTRUCTED SIGNAL ----
    # plt.subplot(313)
    # plt.plot(t, x, label='Original Signal')
    # # Reconstruct signal using interpolated values
    # plt.plot(t, reconstructed, 'r--', label='Reconstructed from Interp. Peak')
    # plt.title('Signal Reconstruction Using Interpolated Peak')
    # plt.xlabel('Time (s)')
    # plt.ylabel('Displacement (m)')
    # plt.legend()

    # # plt.tight_layout()
    # plt.show()
# print(np.round(accuracies*100, 4))
print(f"Average accuracy {round(np.mean(accuracies)*100, 2)}%") 

# RESULT
# excellent accuracies!!! ~99%