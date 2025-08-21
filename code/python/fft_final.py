import numpy as np
from fft_funcs import extracter
import matplotlib.pyplot as plt

# ~~~~~~ Generate fake signal ~~~~~~
sr = 20 #Hz -- sampling rate
# sampling interval
ts = 1.0/sr # timestamp
time_range = 10
t = np.arange(0,time_range,ts)

signals = 1
x = 0
real_power = 0
for i in range(signals):
    frequency = round(2*np.random.rand(), 3)
    amp = round(np.random.rand(), 4)
    x += amp*np.sin(2*np.pi*frequency*t)
    print("Freq %s:" %i, frequency, "Amp %s:" %i, amp)
    real_power += 1/2*((frequency*2*np.pi)**2)*(amp**2)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

num_peaks = 1
amplitudes, frequencies = extracter(x, t, num_peaks) # pass signal, corresponding times, number of peaks to identify

print(f"Detected frequency: {np.round(frequencies, 4)} Hz")
print(f"Peak amplitude: {np.round(amplitudes, 4)} m")


# ---- PLOT ----
plt.plot(t, x, label='Original Signal')
# Reconstruct signal using interpolated values
reconstructed = 0
model = np.zeros(len(t))
modelled_power = 0
for i in range(len(frequencies)):
    reconstructed += amplitudes[i]*np.sin(float(2*np.pi*frequencies[i])*np.array(t))
    modelled_power += 1/2*((frequencies[i]*2*np.pi)**2)*(amplitudes[i]**2)*(4)
plt.plot(t, reconstructed, 'r--', label='Reconstructed from peak amplitude and frequency')
# plt.text(8, 0, f'Power: {round(modelled_power, 2)}W', color='black', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))
plt.title('Signal Reconstruction Using Interpolated Peak')
plt.xlabel('Time (s)')
plt.ylabel('Displacement (m)')
plt.legend(loc='upper right')

# plt.tight_layout()
plt.show()