import numpy as np
import matplotlib.pyplot as plt
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

n = 100
length = 10
t = np.linspace(0, length, n, endpoint=False) # time signal
input_frequency = 2.606
input_amplitude = 1.799
x = input_amplitude*np.sin(input_frequency*2*np.pi* t) + 0.2*np.sin(0.2*2*np.pi* t) # generated signal plus some noise to see how it copes

# ---- PLOT ORGINAL SIGNAL ----
plt.figure(figsize=(12, 8))

plt.subplot(311)
plt.plot(t, x, label='Original Signal')
plt.title('Original Signal')
plt.ylabel('Displacement (m)')
plt.ylim(-4, 4)

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
plt.subplot(312)
plt.stem(freqs, amps, markerfmt=" ", basefmt="-b", label='Zero-padded FFT')
# plt.axvline(freq_peak, color='r', linestyle='--', label=f'Interp. Peak ≈ {freq_peak:.4f} Hz')
plt.title("Frequency Spectrum from FFT")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Amplitude")
plt.legend()

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
for i in range(len(peak_freqs)):
    print(2*np.pi*peak_freqs[i], type(2*np.pi*peak_freqs[i]))
    reconstructed += peak_amps[i]*np.sin(2*np.pi*peak_freqs[i]*t)

# --- PRINT RESULTS ---
print(f"Detected frequency: {np.round(peak_freqs, 4)} Hz")
print(f"Estimated amplitude: {np.round(peak_amps, 4)} m")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ---- PLOT RECONSTRUCTED SIGNAL ----
plt.subplot(313)
plt.plot(t, x, label='Original Signal')
# Reconstruct signal using interpolated values
plt.plot(t, reconstructed, 'r--', label='Reconstructed from Interp. Peak')
plt.title('Signal Reconstruction Using Interpolated Peak')
plt.xlabel('Time (s)')
plt.ylabel('Displacement (m)')
plt.ylim(-4, 4)
plt.legend()

plt.tight_layout()
plt.show()