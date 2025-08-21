# Ocean wave energy device

* [Background](#Background)
* [Power extraction](#Power)
  * [First thoughts and planning](#Intro)
  * [Sensors](#Sensors)
  * [Fast Fourier Transform](#FFT)
  * [Filtering and results](#Filtering)
* [GUI and interface](#GUI)

## Background

This project is to design and test a wave energy converter that will use the motion of waves in the ocean to produce electrical energy. Ocean waves are a huge and largely untapped renewable energy source – wave energy could provide an estimated easily exploitable 500 gigawatts according to [CorPower](http://corpowerocean.com/a-short-history-of-wave-energy/) – and a major advantage of wave energy is that it is not dependent on weather conditions or time of day. Energy sources such as solar power and wind energy rely on variable conditions; solar depends on availability of sunlight, and wind depends on whether and how fast the wind is blowing. In contrast, in many coastal environments waves are constant throughout the day and night, meaning that power can be constantly generated and does not have to be stored to be distributed as much as other renewable energy forms. Wave energy thus has the capacity to provide a stable, reliable source of energy at scale. 

Building on the work done by The Vortical Flow Research Laboratory, this project tests and characterizes an existing wave-energy converter device. The type of wave energy device being studied is a point absorber, which consists of a floating buoy and a fixed rig. The relative motion of the two would drive a power takeoff device (PTO). In our system, the PTO is simulated by a linear damper -- the damper dissipates the energy of the motion of the buoy using electromagnetic resistance. A standard PTO would have a fixed damping coefficient – this setup allows us to test multiple different damping coefficients, simulating different PTO options, without having to reconfigure the device. The goal is to determine the optimal configuration that produces the most power for this particular setup.

<img src="https://github.com/lgray14/Ocean-wave-energy-device/blob/main/images/buoy.jpg" height="300">

<sup> *The device -- a modified Home Depot bucket which has been constrained to move vertically* </sup>


My role in this project has been to run tests, conduct analysis of the data collected, and to add functionality to the device. By adding sensors and writing code, I have been able to calculate and display the instantaneous power output that the device is producing. Through lab tests, I have helped to characterize the linear damper to get an exact relationship between the current to the electromagnets and the damping ratio, and have verified that the code I have worked on functions in lab scenarios. Furthermore, I have integrated my work with the Graphical User Interface (GUI) to better communicate the results. 


## Power

### Intro
To calculate the insantaneous power produced by the device, the motion caused by various signals must be added up. For each signal *j*, the power corresponding to each complex amplitude *X<sub>j</sub>* is given by 
``` math
 P_{j}(\omega) = \frac{1}{2}\omega^{2}b_{PTO}|X_{j}|^{2}
```
So, in order to find the instantaneous power, the complex amplitude (*X<sub>j</sub>* and frequency of the signal (*&omega;*), as well as the damping ratio of the system (*b<sub>PTO</sub>*), must be found. The damping ratio of the system, *b<sub>PTO</sub>*, is, in our case, related to the current being passed to the electromagnets by a constitutive relationship -- Rafael Tejeda has done a bunch of work on characterizing this relationship this summer. Thus it remains to find the amplitude and frequency of the signal. 

Both of these quantities can be found from the position data. A fast Fourier transform (FFT) can be applied to the sinusoidally varying position data to reconstruct from the data points a sinusoidal equation, in particular the amplitude and the frequency. FFTs can also capture the phase but this is not relevant for this data set. 

<img src="https://github.com/lgray14/Ocean-wave-energy-device/blob/main/images/sinusoid.png" height="300">

<sup> *When the bucket moves in the tank, its position varies sinusoidally as shown here on the left hand graph* </sup>


### Sensors

To get the position data, a position sensor must be equipped with the device. The device previously had a distance sensor that worked by turning as the buoy bobbed up and down -- the contact with the frame made it turn and the number of turns corresponded to the distance travelled. However, this sensor added a lot of friction to the system and we decided it was better to use a [time of flight (ToF) sensor](https://www.adafruit.com/product/3317). The VL53L0X ToF sensor seemed like a good choice, since the range need was fairly short and the accuracy need pretty high, and the sensor itself is small and lightweight. The sensor can take readings at different frequencies, so it needed to be tuned to take as many readings as possible in a certain band of accuracy. After some trial, the following Arduino code represents the best balance between speed of readings and accuracy:
```C++
// HIGH ACCURACY VL53L0X CODE

#include <Adafruit_VL53L0X.h>

Adafruit_VL53L0X lox = Adafruit_VL53L0X();

void setup() {
  Serial.begin(115200);

  if (!lox.begin()) {
    Serial.println(F("Failed to boot VL53L0X"));
    while(1);
  }

  // high efficiency setup
  lox.setMeasurementTimingBudgetMicroSeconds(100000);
}

void loop() {
  VL53L0X_RangingMeasurementData_t measure;
  lox.rangingTest(&measure, false); // pass in 'true' to get debug data printout!

  int distMeasured = 0;
  if (measure.RangeStatus != 4) {  // phase failures have incorrect data
    distMeasured = measure.RangeMilliMeter;
  }

  Serial.println(distMeasured);
}
```

The VL530X has been working pretty well, however it is definitely not as accurate as the turning position sensor. However, the friction improvement can be seen very clearly in testing. The friction was decreasing the amplitude of the device's oscillations significantly and was noticeably affecting the damping ratio curve. In particular, the friction was acting as damping of its own and in the free decay tests was creating an undesireable offset in the relationship between current and damping. The ToF sensor is a huge improvement in this area but it definitely has some noticeable noise, in the range of a couple millimeters. It doesn't seem to majorly affect data readings so far, especially given the size of the data sets needed for the FFT process, but this could possibly use some more looking into, and in the future maybe a more precise sensor. 

I placed the sensor on the frame of the wave device to be able to measure the distance to the moving buoy. Having the sensor fixed instead of moving with the system made the wiring less complicated and more secure.

<img src="https://github.com/lgray14/Ocean-wave-energy-device/blob/main/images/wiring1.jpg" height="300"> <img src="https://github.com/lgray14/Ocean-wave-energy-device/blob/main/images/wiring2.jpg" height="300">

The other sensor that I equipped the system with was a Hall Effect sensor, which measures the current being passed to the electromagnets. There is a direct relationship between the current to the electromagnets and the damping ratio because of Lenz's law -- when the magnets are on, eddy currents are generated in the copper plate, leading to a damping effect. It's important to know the current to the magnets for this reason, and the damping ratio goes into the power. The Hall Effect sensor that we chose was the [INA219](https://www.adafruit.com/product/904) because its range of currents was comparable to the testing conditions and for it's lightness and accuracy. 

The INA sensor has also worked pretty well, however it is not quite linear -- it reads slightly too high at low currents and slightly too low at higher currents. However, for our purposes again, this has not been a huge deal but it's not the absolute most accurate. In the long term, having a Hall Effect could be made obsolete by having a tuneable power source. The current DC power supply we have to the electromagnets has to be controlled by turning knobs, but there are fancier to have one that can be controlled from a computer or by a program. Having the ability to tune the current and therefore the damping coefficient could be very useful -- it would make it possible to match specific frequencies with damping coefficients in order to maximize power output.

<img src="https://github.com/lgray14/Ocean-wave-energy-device/blob/main/images/damper.jpg" height="300">

<sup> *Damper system* </sup>

The following is the Arduino code for the ToF and Hall Effect sensors, as well as the force probe the system was already equipped with.

```C++
#include "HX711.h"
#include <Wire.h>
#include <Adafruit_INA219.h>
#include <Adafruit_VL53L0X.h>
#include <Simpletimer.h>

//hall effect sensor
Adafruit_INA219 ina219;
// ToF sensor
Adafruit_VL53L0X lox = Adafruit_VL53L0X();
// timer
Simpletimer timer1{};
// Load Cell Amp -- force readings
#define calibration_factor 1250.0 //  This value is obtained using the SparkFun_HX711_Calibration sketch
#define DOUT  3
#define CLK  2
HX711 scale;


void setup() {
  // initialize serial communication at 9600 bits per second:
  Serial.begin(115200);

  // initialize load cell amp
  scale.begin(DOUT, CLK);
  scale.set_scale(calibration_factor); //This value is obtained by using the SparkFun_HX711_Calibration sketch
  scale.tare(); //Assuming there is no weight on the scale at start up, reset the scale to 0
  
  // initialize Hall Effect sensor
  ina219.begin();
  
  // set up ToF sensor
  if (!lox.begin()) {
    Serial.println(F("Failed to boot VL53L0X"));
    while(1);
  }
  lox.setMeasurementTimingBudgetMicroSeconds(100000); // from testing, this value results in a good balance between noise and sampling rate
}

// set sampling rate -- should not change
float sampling_rate = 10; //Hz
float gap = (1/sampling_rate)*1000;

void loop() {
  // ToF sensor
  VL53L0X_RangingMeasurementData_t measure;
  lox.rangingTest(&measure, false); // pass in 'true' to get debug data printout!

  if (timer1.timer(gap)) { // take samples at specified rate
    int distMeasured = 0;
    if (measure.RangeStatus != 4) {  // phase failures have incorrect data
      distMeasured = measure.RangeMilliMeter;
    }
    // Hall effect sensor
    float current_mA = 0;
    current_mA = ina219.getCurrent_mA();

    Serial.print(distMeasured);
    Serial.print(' ');
    // SPACER FOR GUI READINGS
    Serial.print(0);
    Serial.print(' ');
    // 
    Serial.print(current_mA, 2); // current returns a float ==> GUI does not cope well, resolve
    Serial.print(' ');
    Serial.println(scale.get_units(), 2); //force scale.get_units() returns a float
  }
}
```

### FFT

With the sensors attached and the position in particular readable from the Arduino, I turned my attention to extracting the amplitude and frequency from the signal. My first step was actually a pretty signficant misstep -- I assumed the amplitude would be pretty variable so I opted for a version of the fast Fourier transform (FFT) called a Hilbert-Huang transform. The Hilbert transform is great at extracting the instantaneous amplitude and amplitude envelope, and at being very accurate with a smaller data set, but lacks in the frequency extraction. It's much better at finding the natural frequency of the decaying system, and doesn't do a great job with a consistent forcing frequency. After realizing this, I switched over to a traditional FFT model. FFT does a great job with a consistent frequency and amplitude, which is ideal for the consistent signals we are working with.

The fast Fourier transform takes a set of data points and extracts the amplitude and frequency of the sinusoid it can recognize in that data. For this reason it works best with a couple oscillations. By trying a couple different time intervals/data set lengths, it seems like at the low frequencies that the wave tank usually entails that about ten seconds worth of data, at around 100 points, is ideal. The next problem that comes up from this is how to have a moving data set of 100 points. The power should be live and change with the signal, so the last ten seconds worth of data should be the ones being evaluated. By storing the data and taking the last 100 points, this can be achieved.

The next issue I dealt with was low sampling frequency. The ToF sensor has a pretty high sampling frequency but unless I made the FFT time interval much longer, the resolution for the detectable frequencies was too low. With the setup I had, the most accuracy I could get was .1 Hz, which isn't terrible, but even small variations in the frequency of the incoming waves make large differences in the output of the device. In order to resolve this, I used zero-padding and windowing. I used Hann windowing specifically, which smooths the data and helps the algorithm find the peaks that the data might not quite show. At first, I tried various types of interpolation, and then zero-padding on its own, but I got by far the best results by using windowing and zero-padding together. These are included in the "extracter" function shown below. After that, I also used a quadratic interpolation to further smooth and identiify the peaks. This is potentially a little overkill but the results are very consistent and it doesn't take an enormous amount of time. 

**FFT enabling functions**
```python
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
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

def extracter(x, t, k):
    """
    Extract k peak ampliudes and their corresponding frequencies using zero-padded fast Fourier transform.

    Parameters:
        x: original displacement signal
        t: corresponding time data
        k: number of peak frequencies to be identified

    Returns:
        peak_amps: identified peak amplitudes, list of length k
        peak_freqs: corresponding frequencies
    """

    n = len(t)
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

    # find peaks to analyze around
    peaks, _ = find_peaks(amps)
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
    
    return peak_amps, peak_freqs
```

### Filtering

Lastly, the final major hurdle was dealing with layered signals. For the purposes of the wave tank, identifying the most dominant signal is the priority, but I wanted to include functionality to expand a little bit how many peaks can be identified. The "extracter" function is arranged to take *k* as an input -- the function will return the *k* largest peaks. This is shown in the image below on the left, where the data contains two signals but the result filters to the most dominant one. Since this was not a top priority it's not super robust, so you do have to identify a number of peaks to seek out by entering a value *k*. Additionally, the windowing method does lead to 'echo peaks' (not sure if there is a technical term for this) where around the actual peak there are several short rebound peaks that can pop up. This effect can be seen in the image below on the right for a two-signal data set -- there are two real peaks, but around the base of the taller one there are a couple little bumps. As of yet I have not figured out how to get the peak finder to ignore these as they cannot be distinguished from the real peaks except by their proximity to one and that they are symmetric about a larger one. Again this was a lower priority and can be improved in the future. Another future project that would improve the effectiveness of the program would be to auomatically identify all the major peaks instead of having to enter a value *k* for how many peaks to look for.

<img src="https://github.com/lgray14/Ocean-wave-energy-device/blob/main/images/2signal_fft.png" height="300"> <img src="https://github.com/lgray14/Ocean-wave-energy-device/blob/main/images/fft_echos.png" height="300">

The following is an example of how the functions can be implemented to analyze a fake signal. 

**Final FFT example**
```python
import numpy as np
from fft_funcs import extracter
import matplotlib.pyplot as plt

# ~~~~~~ Generate fake signal ~~~~~~
sr = 10 #Hz -- sampling rate
# sampling interval
ts = 1.0/sr # timestamp
time_range = 10
t = np.arange(0,time_range,ts)

signals = 2
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
plt.text(4.5, 0, f'Power: {round(modelled_power, 2)}W', color='black', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))
plt.title('Signal Reconstruction Using Interpolated Peak')
plt.xlabel('Time (s)')
plt.ylabel('Displacement (m)')
plt.legend()

# plt.tight_layout()
plt.show()
```

**Result of this code with two signals**

<img src="https://github.com/lgray14/Ocean-wave-energy-device/blob/main/images/2signal.png" height="500">

## GUI
