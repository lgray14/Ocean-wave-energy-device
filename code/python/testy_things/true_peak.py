import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
from numpy.random import rand
from scipy.interpolate import interp1d
from scipy.stats import zscore

sr = 60 #Hz -- sampling rate
# sampling interval
ts = 1.0/sr # timestamp
time_range = 10
t = np.arange(0,time_range,ts)

faccuracy = []
aaccuracy = []
paccuracy = []
for j in range(1000):
    print(j)
    # Generate fake signal
    n = 2
    x = 0
    real_power = 0
    input_amps = []
    input_freqs = []
    for i in range(n):
        frequency = round(sr/4*rand(), 2)
        amp = round(rand(), 2)
        input_freqs.append(frequency)
        input_amps.append(amp)
        x += amp*np.sin(2*np.pi*frequency*t)
        # print("Freq %s:" %i, frequency, "Amp %s:" %i, amp)
        real_power += 1/2*((frequency*2*np.pi)**2)*(amp**2)
    # plt.figure(figsize = (9, 5))
    # plt.subplot(311)
    # plt.plot(t, x, 'r')
    # plt.ylabel('Amplitude (m)')
    # plt.xlabel('Time(s)')

    X=np.fft.fft(x)
    freq = np.fft.fftfreq(len(t), 1/sr)

    n_oneside = int(len(freq)//2)
    f_oneside = freq[:n_oneside]

    X_oneside = abs(X[:n_oneside]/n_oneside)

    # plt.subplot(312)
    # plt.stem(f_oneside, X_oneside, 'b', markerfmt=" ", basefmt="-b")
    # plt.xlabel('Freq (Hz)')
    # plt.ylabel('Normalized FFT Amplitude')

    peaks, _ = find_peaks(X_oneside, threshold=.005)
    # dom_amps = X_oneside[peaks]
    # dom_freqs = f_oneside[peaks]

    true_amps = []
    true_freqs = []
    if len(peaks) > 0 and min(f_oneside[peaks]) > .2: # throws errors if no peaks to work with or peak freq is too small
        for i in range(len(peaks)):
            f_peak_region = f_oneside[peaks[i]-3:peaks[i]+4]  # Select a small region around the peak
            X_peak_region = X_oneside[peaks[i]-3:peaks[i]+4]

            # Interpolate the region (e.g., using cubic interpolation)
            interp_func = interp1d(f_peak_region, X_peak_region, kind='cubic')
            f_fine = np.linspace(f_peak_region[0], f_peak_region[-1], 1000)
            X_fine = interp_func(f_fine)

            # plt.plot(f_fine, X_fine, 'b--')

            # Find the refined peak
            refined_peak_index = np.argmax(X_fine)
            true_freqs.append(f_fine[refined_peak_index])
            true_amps.append(X_fine[refined_peak_index])
        # print("Identified peak amplitudes:", true_amps)
        # print("Identified corresponding frequencies:", true_freqs)

        # # plt.subplot(313)
        # model = np.zeros(len(t))
        modelled_power = 0
        for i in range(len(true_amps)):
            modelled_power += 1/2*((true_freqs[i]*2*np.pi)**2)*(true_amps[i]**2)
            # for j in range(len(t)):
            #     model[j] += true_amps[i]*np.sin(2*np.pi*true_freqs[i]*t[j])
        # plt.plot(t, model, 'g')
        # plt.ylabel("Modelled Amplitude (m)")
        # plt.xlabel("Time (s)")

        # print("Real Power up to damping coefficient:", real_power)
        # print("Modelled Power up to damping coefficient:", modelled_power)
        if real_power > 0:
            if modelled_power > 0:
                faccuracy.append(max(true_freqs)/max(input_freqs))
                aaccuracy.append(max(true_amps)/max(input_amps))
                paccuracy.append(modelled_power/real_power)
        # plt.show()

# f_zscore = zscore(faccuracy)
# a_zscore = zscore(aaccuracy)
# p_zscore = zscore(paccuracy)
# filt_facc = [x for x, z in zip(faccuracy, f_zscore) if abs(z) < 2]
# filt_aacc = [x for x, z in zip(aaccuracy, a_zscore) if abs(z) < 2]
# filt_pacc = [x for x, z in zip(paccuracy, p_zscore) if abs(z) < 2]

# # print([round(float(num), 2) for num in filt_facc])
# print("Average frequency accuracy:", np.mean(filt_facc))
# # print([round(float(num), 2) for num in filt_aacc])
# print("Average amplitude accuracy:", np.mean(filt_aacc))
# # print([round(float(num), 2) for num in filt_pacc])
# print("Average power accuracy:", np.mean(filt_pacc))

# unfiltered
print("Average frequency accuracy:", np.mean(faccuracy))
print("Average ampltiude accuracy:", np.mean(aaccuracy))
print("Average power accuracy:", np.mean(paccuracy))

# QUADRATIC RESULTS: 
# average accuracy of freqquency: ~ 95%
# average accuracy of amplitude ~ 87% 
# average accuracy of power: ~76%

# CUBIC RESULTS:
# Average frequency accuracy: ~95%
# Average amplitude accuracy: ~87%
# Average power accuracy: ~76%

# SLINEAR
# Average frequency accuracy: ~95%
# Average amplitude accuracy: ~86%
# Average power accuracy: ~76%

# exact same for everything -- what on earth
# also exactly the same as with no interp :((((((((