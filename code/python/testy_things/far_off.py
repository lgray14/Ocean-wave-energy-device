## TEST HOW FAR OFF POWER IS FROM ACTUAL

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
from numpy.random import rand

percent_accuracy = np.zeros(1000)

for i in range(1000):
    sr = 60 #Hz -- sampling rate
    # sampling interval
    ts = 1.0/sr # timestamp
    time_range = 10
    t = np.arange(0,time_range,ts)

    frequency = round(sr/8*rand(), 2)
    amp = round(rand(), 2)
    x = amp*np.sin(2*np.pi*frequency*t)
    real_power = 1/2*((frequency*2*np.pi)**2)*(amp**2)
    # print(real_power)

    # frequency = round(sr/8*rand(), 2)
    # amp = round(rand(), 2)
    # x += amp*np.sin(2*np.pi*frequency*t)
    # real_power += 1/2*((frequency*2*np.pi)**2)*(amp**2)

    # frequency = round(sr/4*rand(), 2)
    # amp = round(rand(), 2)
    # x += amp*np.sin(2*np.pi*frequency*t)
    # real_power += 1/2*((frequency*2*np.pi)**2)*(amp**2)

    # frequency = round(sr/4*rand(), 2)
    # amp = round(rand(), 2)
    # x += amp*np.sin(2*np.pi*frequency*t)
    # real_power += 1/2*((frequency*2*np.pi)**2)*(amp**2)

    # do FFT
    X=np.fft.fft(x)
    freq = np.fft.fftfreq(len(t), 1/sr)
    # print(freq)

    # peaks can only be validly found at half the sampling frequency -- others reflect: only process the first half of frequencies
    n_oneside = int(len(freq)//2)
    # get the one side frequency
    f_oneside = freq[:n_oneside] # limit to the truly observable frequencies -- half of the sampling frequency

    # normalize the amplitude
    X_oneside = abs(X[:n_oneside]/n_oneside)
    # print(X_oneside, len(X_oneside))

    peaks, _ = find_peaks(X_oneside, threshold=.005)
    # print(peaks)
    # print("Amplitudes:", X_oneside)
    dom_amps = X_oneside[peaks]
    dom_freqs = f_oneside[peaks]
    # print("Dominant frequencies:", dom_freqs, "   Corresponding amplitudes:", dom_amps)

    # plot modelled curve
    model = np.zeros(len(t))
    modelled_power = 0
    for j in range(len(dom_amps)):
        modelled_power += 1/2*((dom_freqs[j]*2*np.pi)**2)*(dom_amps[j]**2)
    
    print("\nTrial:", i)
    print(real_power)
    print(modelled_power)

    if real_power > 0:
        if modelled_power > 0:
            percent_accuracy[i] = modelled_power/real_power
        else:
            pass
    # print("Trial:", i, "Accuracy:", percent_accuracy)

print(percent_accuracy)
mean_accuracy = np.mean(percent_accuracy)
print("Average accuracy of power calculation vs real", mean_accuracy)

## RESULT: POWER IS AN AVERAGE OF ~75% OF ACTUAL :(
#  ==> Case where resolution is .1 Hz and input signal resolution is .01 Hz 
#       -- should ask about how accurate is wanted