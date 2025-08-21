import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
from numpy.random import rand
from scipy.interpolate import interp1d
from scipy.stats import gaussian_kde

sr = 20 #Hz -- sampling rate
# sampling interval
ts = 1.0/sr # timestamp
time_range = 10
t = np.arange(0,time_range,ts)

# Generate fake signal
n = 1
x = 0
real_power = 0
for i in range(n):
    frequency = round(sr/4*rand(), 2)
    amp = round(rand(), 2)
    x += amp*np.sin(2*np.pi*frequency*t)
    print("Freq %s:" %i, frequency, "Amp %s:" %i, amp)
    real_power += 1/2*((frequency*2*np.pi)**2)*(amp**2)
plt.figure(figsize = (9, 5))
plt.subplot(311)
plt.plot(t, x, 'r')
plt.ylabel('Amplitude (m)')
plt.xlabel('Time(s)')

X=np.fft.fft(x)
freq = np.fft.fftfreq(len(t), 1/sr)

n_oneside = int(len(freq)//2)
f_oneside = freq[:n_oneside]

X_oneside = abs(X[:n_oneside]/n_oneside)

plt.subplot(312)
plt.stem(f_oneside, X_oneside, 'b', markerfmt=" ", basefmt="-b")
plt.xlabel('Freq (Hz)')
plt.ylabel('Normalized FFT Amplitude')

peaks, _ = find_peaks(X_oneside, threshold=.005)
dom_amps = X_oneside[peaks]
dom_freqs = f_oneside[peaks]
print("Dominant frequencies:", dom_freqs, "   Corresponding amplitudes:", dom_amps)

# true_amps = []
# true_freqs = []
# for i in range(len(peaks)):
#     f_peak_region = f_oneside[peaks[i]-7:peaks[i]+8]  # Select a small region around the peak
#     X_peak_region = X_oneside[peaks[i]-7:peaks[i]+8]

#     # Interpolate the region (e.g., using cubic interpolation)
#     interp_func = interp1d(f_peak_region, X_peak_region, kind='cubic')
#     f_fine = np.linspace(f_peak_region[0], f_peak_region[-1], 1000)
#     X_fine = interp_func(f_fine)

#     plt.plot(f_fine, X_fine, 'b--')

#     # Find the refined peak
#     refined_peak_index = np.argmax(X_fine)
#     true_freqs.append(f_fine[refined_peak_index])
#     true_amps.append(X_fine[refined_peak_index])

    # plt.plot(f_vals, pdf_vals, "r--")

def parabolic_interpolation(freqs, magnitudes):
    idx = np.argmax(magnitudes)
    if idx <= 0 or idx >= len(magnitudes) - 1:
        return freqs[idx], magnitudes[idx]
    
    y0, y1, y2 = magnitudes[idx-1], magnitudes[idx], magnitudes[idx+1]
    x0, x1, x2 = freqs[idx-1], freqs[idx], freqs[idx+1]
    delta = 0.5 * (y0 - y2) / (y0 - 2*y1 + y2)
    freq_peak = x1 + delta * (x2 - x1)
    mag_peak = y1 - 0.25 * (y0 - y2) * delta
    return freq_peak, mag_peak
true_amps, true_freqs = parabolic_interpolation(f_oneside, X_oneside)
print(true_amps, true_freqs)

# print("Identified peak amplitudes:", [float(num) for num in true_amps], "  Identified corresponding frequencies:", [float(num) for num in true_freqs])
print("Identified peak amplitudes:", true_amps, "  Identified corresponding frequencies:", true_freqs)

plt.subplot(313)
model = np.zeros(len(t))
modelled_power = 0
for i in range(len(true_amps)):
    modelled_power += 1/2*((true_freqs[i]*2*np.pi)**2)*(true_amps[i]**2)
    for j in range(len(t)):
        model[j] += true_amps[i]*np.sin(2*np.pi*true_freqs[i]*t[j])
plt.plot(t, model, 'g')
plt.ylabel("Modelled Amplitude (m)")
plt.xlabel("Time (s)")

print("Real Power up to damping coefficient:", real_power)
print("Modelled Power up to damping coefficient:", modelled_power)
print("Power accuracy: %s percent" % ((modelled_power/real_power)*100))

plt.show()