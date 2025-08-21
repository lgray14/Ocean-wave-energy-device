import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import fft, fftfreq
from scipy.signal import find_peaks


def sinc_interp(x, s, u):
    """
    Perform sinc interpolation.

    Parameters:
        x (ndarray): Original signal values (samples).
        s (ndarray): Sample points corresponding to `x`.
        u (ndarray): Desired interpolation points.

    Returns:
        ndarray: Interpolated signal at points `u`.
    """
    # Compute the sinc matrix
    sinc_matrix = np.sinc((u[:, None] - np.array(s)[None, :])/(s[1] - s[0]))
    # Perform interpolation
    results = np.dot(sinc_matrix, x)
    # print(results)
    # interp_mag = np.abs(results)/(n//2) # normalize interpolation
    return results

# def sinc_interp(x, s, u):
#     """Sinc interpolation of signal s (sampled at x) evaluated at points u."""
#     # Broadcasting s over u-x. Use outer difference.
#     # print(s.shape)
#     # print("u:", u, u.shape, np.array(u)[:, None])
#     # print("x:", x, np.array(x).shape, np.array(x)[None, :])
#     # print(np.sinc((np.array(u)[:, None] - np.array(x)[None, :]) / (x[1]-x[0])))
    
#     # return np.dot(s, np.sinc((u[:, None] - x[None, :]) / (x[1]-x[0])))
#     return np.dot(s, np.sinc((np.array(u)[:, None] - np.array(x)[None, :]) / (x[1]-x[0])))

# Parameters
n = 500 # number of points
length = 10
# omega = 2.0 * np.pi / Lx
x = np.linspace(0, length, int(n), endpoint=False)

# Composite signal: two sinusoids, one with fractional freq
y = 2.0 * np.sin(2.43 * 2*np.pi * x) + .2 * np.sin(.2 * 2*np.pi * x)
plt.subplot(411)
plt.plot(x,y)
plt.xlabel('Time(s)')
plt.ylabel('Displacement (m)')
plt.ylim(-4, 4)

# FFT
fft_vals = fft(y)
freqs_raw = fftfreq(int(n), d=x[1]-x[0])
# print(freqs_raw, len(freqs_raw))
mask = freqs_raw >= 0
freqs = freqs_raw[mask]
freqs = [float(freq) for freq in freqs]
# print(freqs, len(freqs))
fft_vals_half = fft_vals[:n//2]
amps = np.abs(fft_vals[:n//2]/(n//2))  # single-sided amplitude

peaks, _ = find_peaks(amps, threshold=.005)
# print(freqs, len(freqs))
# print(amps, len(amps))
print(peaks)
# print("Amplitudes:", X_oneside)
dom_amps = amps[peaks]
# print("Dominant amps:", dom_amps)
dom_freqs = np.array(freqs)[peaks]
print("PRELIMINARY")
print("Dominant frequencies:", dom_freqs, "   Corresponding amplitudes:", dom_amps)

# filter out k highest peaks

k = 1
combined_lists = list(zip(amps, freqs))
sorted_lists = sorted(combined_lists, key=lambda x: x[0]) # sort BOTH frequencies and amplitudes based on amplitude

filter_amps = [item[0] for item in sorted_lists][-k:]
# print(sorted(amps), filter_amps)
filter_freqs = [item[1] for item in sorted_lists][-k:]
# print(sorted(freqs), filter_freqs)

plt.subplot(412)
plt.plot(freqs, amps, 'o', markersize=3, label='Original FFT bins')
plt.xlabel('Frequency')
plt.ylabel('Amplitude')
plt.legend()
plt.grid(True)
# interpolate around peak
for i in range(k):
    dense_freqs = np.linspace(filter_freqs[i]-.25, filter_freqs[i]+.25, 1000)
    # print(dense_freqs)
    # print(len(freqs))
    # print(len(amps))
    # print(len(dense_freqs))
    interp_fft = sinc_interp(np.log(amps), freqs, dense_freqs) # interpolate complex fft vals ==> output is normalized mags
# Find peak around 5.4 Hz
# interp_mag = np.abs(interp_fft)/(n//2)
interp_mag = interp_fft

peak_idx = np.argmax(interp_mag)
freq_peak = dense_freqs[peak_idx]
amp_peak = interp_mag[peak_idx]
phase_peak = np.angle(interp_fft[peak_idx])

print(f"Detected peak at ~{freq_peak:.4f} Hz with amplitude ~{amp_peak:.3f}")

# Plot
# plt.figure(figsize=(8,4))
# plt.plot(dense_freqs, interp_mag, label='Interpolated FFT')

plt.subplot(413)
plt.title('Reconstructed Signal')
reconstructed = amp_peak * np.sin(2*np.pi*freq_peak*x + phase_peak)
plt.plot(x, reconstructed, 'r', label='Sinc (amp + phase)')
plt.plot(x, max(filter_amps)*np.sin(2*np.pi*max(filter_freqs)*x), label='Simple peak')
plt.xlabel('Time (s)')
plt.ylabel('Displacement (m)')
plt.ylim(-4, 4)
plt.legend()

# plt.subplot(414)
# plt.plot(x, y) # original
# # plt.plot(x, (interp_mag[peak_idx]*np.sin(2*np.pi*dense_freqs[peak_idx]*x)), 'r', label='Sinc')
# plt.plot(x, (max(filter_amps)*np.sin(2*np.pi*max(filter_freqs)*x)), label='Simple peak')
# # plt.plot(x, (max(filter_amps)*np.sin(2*np.pi*dense_freqs[peak_idx]*x)), 'g', label='combined')

plt.show()