import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

# # # Generate a sample signal: a combination of two sine waves
# sampling_rate = 1000  # Hz
# duration = 1.0  # seconds
# t = np.linspace(0, duration, int(sampling_rate * duration), endpoint=False)

# # Frequencies of the sine waves
# freq1 = 50  # Hz
# freq2 = 120  # Hz

# # Create the signal
# signal = np.sin(2 * np.pi * freq1 * t) + 0.5 * np.sin(2 * np.pi * freq2 * t)

# # Perform FFT
# fft_result = np.fft.fft(signal)
# frequencies = np.fft.fftfreq(len(t), 1 / sampling_rate)

# # Take the magnitude of the FFT result
# magnitude = np.abs(fft_result)

# # Plot the original signal
# plt.figure(figsize=(12, 6))
# plt.subplot(2, 1, 1)
# plt.plot(t, signal)
# plt.title("Original Signal")
# plt.xlabel("Time (s)")
# plt.ylabel("Amplitude")

# # Plot the FFT result (only positive frequencies)
# plt.subplot(2, 1, 2)
# plt.plot(frequencies[:len(frequencies)//2], magnitude[:len(magnitude)//2])
# plt.title("FFT of the Signal")
# plt.xlabel("Frequency (Hz)")
# plt.ylabel("Magnitude")
# plt.tight_layout()
# plt.show()

def sinc_interp(x, s, t):
    """
    Perform 1D sinc interpolation.
    
    Parameters:
        x : array-like
            The sampled signal values.
        s : array-like
            The sample points corresponding to `x`.
        t : array-like
            The desired interpolation points.
    
    Returns:
        y : array-like
            Interpolated signal values at `t`.
    """
    sinc_matrix = np.sinc((t[:, None] - s[None, :]) / (s[1] - s[0]))
    return np.dot(sinc_matrix, x)

# Original sampled signal
s = np.linspace(0, 10, 11)  # Sample points
# print(s)
x = np.sin(s)               # Sampled signal values
# print(x)

# Interpolation points
interp = np.linspace(0, 10, 100)  # Desired points for interpolation
# print(t)

# Perform sinc interpolation
y = sinc_interp(x, s, interp)
# peaks, _ = find_peaks(y, threshold=.01)
print(y)
# print(y[peaks])

# Plot the results
plt.figure(figsize=(8, 4))
plt.plot(s, x, 'o', label='Sampled Points', color='red')
# plt.plot(s[peaks], y[peaks], 'x', label='Peaks', color = 'green')
plt.plot(interp, y, '-', label='Sinc Interpolation', color='blue')
plt.legend()
plt.title("1D Sinc Interpolation")
plt.xlabel("t")
plt.ylabel("Signal")
plt.grid()
plt.show()