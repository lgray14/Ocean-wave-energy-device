import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from numpy.random import rand, randint

plt.figure(figsize = (10, 6))
plt.subplot(3, 1, 1)
plt.xlabel('Time(s)')
plt.ylabel('Heave distance (m)')

thresh_right = 0
prom_right = 0

for i in range(100):
    t = []
    z = []

    sr = 128 #Hz -- sampling rate
    # sampling interval
    ts = 1.0/sr # timestamp
    time_range = 2
    t = np.arange(0,time_range,ts) # time points 0-4s

    freq = round(sr/4*rand(), 2)
    amp = round(rand(), 2)
    z = amp*np.sin(2*np.pi*freq*t)

    freq = round(sr/4*rand(), 2)
    amp = round(rand(), 2)
    z += amp*np.sin(2*np.pi*freq*t)

    freq = round(sr/2*rand(), 2)
    amp = round(rand(), 2)
    z += amp*np.sin(2*np.pi*freq*t)

    freq = round(sr/2*rand(), 2)
    amp = round(rand(), 2)
    z += amp*np.sin(2*np.pi*freq*t)

    Z=np.fft.fft(z)

    # calculate the frequency
    N = len(Z)
    # print(N)
    n = np.arange(N)
    # print(n)
    T = N/sr
    freq = n/T 
    # print(freq, len(freq))

    # peaks can only be validly found at half the sampling frequency -- others reflect: only process the first half of frequencies
    n_oneside = int(len(freq)/2)
    # print("n_oneside:", n_oneside)
    # get the one sided frequency
    f_oneside = freq[:n_oneside] # limit to the truly observable frequencies -- half of the sampling frequency
    # print("Max frequency:", f_oneside[-1], f_oneside)

    # normalize the amplitude
    Z_oneside = abs(Z[:n_oneside]/n_oneside)

    thresh_peaks, _ = find_peaks(Z_oneside, threshold=.005)
    prom_peaks, _ = find_peaks(Z_oneside, prominence=.01)
    # print(peaks)
    # print("Amplitudes:", X_oneside)
    dom_amps = Z_oneside[thresh_peaks]
    dom_freqs = f_oneside[thresh_peaks]
    if len(dom_amps) == 4:
        thresh_right += 1
    print("\nPEAKS FROM THRESHOLD")
    print("Dominant frequencies:", dom_freqs, "   Corresponding amplitudes:", dom_amps)
    dom_amps = Z_oneside[prom_peaks]
    dom_freqs = f_oneside[prom_peaks]
    if len(dom_amps) == 4:
        prom_right += 1
    print("PEAKS FROM PROMINENCE")
    print("Dominant frequencies:", dom_freqs, "   Corresponding amplitudes:", dom_amps)

print("Prominence accuracy: %s percent" % prom_right)
print("Threshold accuracy: %s percent" % thresh_right)

plt.plot(t, z, 'r')

plt.subplot(3, 1, 2)
plt.stem(f_oneside, Z_oneside, 'b', markerfmt=" ", basefmt="-b")
plt.xlabel('Freq (Hz)')
plt.ylabel('Normalized FFT Amplitude')
# plt.tight_layout()

# plot modelled curve
plt.subplot(3, 1, 3)
model = np.zeros(len(t))
for i in range(len(dom_amps)):
    for j in range(len(t)):
        model[j] += dom_amps[i]*np.sin(2*np.pi*dom_freqs[i]*t[j])
plt.plot(t, model, 'g')
plt.ylabel("Modelled Amplitude (m)")
plt.xlabel("Time (s)")

plt.show()

# results: low prominence does not in fact throw errors readily. it's comparable to/slightly better than threshold -- 
# though very low threshold also does not really throw errors. end result -- pretty much comparable/similar accuracy for 
# use case. 