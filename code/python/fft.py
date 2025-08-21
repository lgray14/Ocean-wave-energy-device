import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
from numpy.random import rand

sr = 50 #Hz -- sampling rate
# sampling interval
ts = 1.0/sr # timestamp
time_range = 10
t = np.arange(0,time_range,ts)

# Generate fake signal
n = 2
x = 0
real_power = 0
for i in range(n):
    frequency = round(sr/4*rand(), 2)
    amp = round(rand(), 2)
    x += amp*np.sin(2*np.pi*frequency*t)
    print("Freq %s:" %i, frequency, "Amp %s:" %i, amp)
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

plt.figure(figsize = (9, 5))
plt.subplot(311)
plt.plot(t, x, 'r')
plt.ylabel('Amplitude (m)')
plt.xlabel('Time(s)')

# do FFT
X=np.fft.fft(x)
# FED FREQUENCY
freq = np.fft.fftfreq(len(t), 1/sr)
# GUI RESULTS
# freq = np.fft.fftfreq(len(t), avg_gap)
# print(freq)

# peaks can only be validly found at half the sampling frequency -- others reflect: only process the first half of frequencies
n_oneside = int(len(freq)//2)
# get the one side frequency
f_oneside = freq[:n_oneside] # limit to the truly observable frequencies -- half of the sampling frequency

# normalize the amplitude
X_oneside = abs(X[:n_oneside]/n_oneside)
# print(X_oneside, len(X_oneside))

plt.subplot(312)
plt.stem(f_oneside, X_oneside, 'b', markerfmt=" ", basefmt="-b")
plt.xlabel('Freq (Hz)')
plt.ylabel('Normalized FFT Amplitude')
plt.tight_layout()

peaks, _ = find_peaks(X_oneside, threshold=.005)
# print(peaks)
# print("Amplitudes:", X_oneside)
dom_amps = X_oneside[peaks]
dom_freqs = f_oneside[peaks]
print("PRELIMINARY")
print("Dominant frequencies:", dom_freqs, "   Corresponding amplitudes:", dom_amps)

# filter out k highest peaks
k = 1
combined_lists = list(zip(dom_amps, dom_freqs))
sorted_lists = sorted(combined_lists, key=lambda x: x[0]) # sort BOTH frequencies and amplitudes based on amplitude

filter_amps = [item[0] for item in sorted_lists][-k:]
# print(sorted(dom_amps), filter_amps)
filter_freqs = [item[1] for item in sorted_lists][-k:]
# print(sorted(dom_freqs), filter_freqs)

# plot modelled curve
plt.subplot(313)
model = np.zeros(len(t))
modelled_power = 0
for i in range(len(filter_amps)):
    modelled_power += 1/2*((filter_freqs[i]*2*np.pi)**2)*(filter_amps[i]**2)
    for j in range(len(t)):
        model[j] += filter_amps[i]*np.sin(2*np.pi*filter_freqs[i]*t[j])
plt.plot(t, model, 'g')
plt.ylabel("Modelled Amplitude (m)")
plt.xlabel("Time (s)")

# print("Real Power up to damping coefficient:", real_power)
print("Modelled Power up to damping coefficient:", modelled_power, "m^2*s^2")

plt.show()