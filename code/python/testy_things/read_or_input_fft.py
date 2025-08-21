import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from numpy.random import rand, randint

t = []
z = []

plt.figure(figsize = (10, 6))
plt.subplot(3, 1, 1)
plt.xlabel('Time(s)')
plt.ylabel('Heave distance (m)')

#read data from GUI fake reader
#--------------------------------------------------------
# with open("gui_results_test.txt", "r") as f:
#     lines = f.readlines()

#     # Loop through all lines, ignoring header.
#     # Add last element to list (i.e. the process name)
#     for l in lines[2:]:
#         z.append(float(l.split()[1]))
#         t.append(float(l.split()[0]))
# plt.plot(t, [item/100 for item in z], 'r') # convert to m rather than cm

# # find average time gap ==> sampling frequency
# diffs = []
# for i in range(len(t)-1):
#     diffs.append(t[i+1] - t[i])
# avg_gap = np.mean(diffs)
# # print(avg_gap)
# sr = 1/avg_gap
#--------------------------------------------------------

# print("Time [s]", t)
# print("Heave distance [cm]", z)


# direct input
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

plt.plot(t, z, 'r')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

plt.subplot(3, 1, 2)
plt.stem(f_oneside, Z_oneside, 'b', markerfmt=" ", basefmt="-b")
plt.xlabel('Freq (Hz)')
plt.ylabel('Normalized FFT Amplitude')
# plt.tight_layout()


# different filtering conditions:
# height: only counts as a peak if it is a certain height 
#  ==> eg, amplitude of .1m: misses lower peaks and is very rigid, doesn't adapt with short waves
# prominence: only counts as a peak if it stands out from data 
#  ==> eg, peak has prominence of .25: theoretically good for noise filtration but might cut out peaks if they are close together
# threshold: only counts a peak if it is a certain height relatve to surrounding data
#  ==> eg, peak is .05m above data point to left and right: might be a good balance but could also miss low peaks or hear big noise

peaks, _ = find_peaks(Z_oneside, threshold=.05)
# print(peaks)
# print("Amplitudes:", X_oneside)
dom_amps = Z_oneside[peaks]
dom_freqs = f_oneside[peaks]
print("Dominant frequencies:", dom_freqs, "   Corresponding amplitudes:", dom_amps)

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