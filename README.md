# Ocean wave energy device

* [Background](#Background)
* [Power extraction](#Power)
  * [First thoughts and planning](#Intro)
  * [Sensors](#Sensors)
  * [Fast Fourier Transform](#FFT)
  * [Filtering and results](#Filtering)
* [GUI and interface](#GUI)
* [Final thoughts](#Final)

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

The following is an example of how the functions can be implemented to analyze a live signal. 

**Final FFT example**
```python
import serial
from time import time
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
from numpy.fft import fft, fftfreq
from scipy.signal import windows, find_peaks
from fft_funcs import extracter

ser = serial.Serial('COM7', baudrate=115200)
print(ser.name)

x = []
t = []
I = []

num_peaks = 2
times = 0

ser.readline() # read and do nothing with first line
start = time()
while True:
    line = ser.readline().decode('utf-8').strip()  # Read a line and decode it
    print(line)
    if not line:  # Stop if no data is received
        break

    x.append(float(line.split()[0])/1000) # convert mm to m
    t.append((time() - start))
    I.append(float(line.split()[2])/1000) # convert from mA to A
    # print(x)
    # print(len(x))
    # print(len(x)%10)
    if len(x) >= 100 and len(x)%10 == 0:
        x100 = x[-100:]
        t100 = t[-100:]
        current = np.mean(I[-100:])
        
        # print(np.round(x100, 3), len(x100))
        # print(np.round(t100, 3), len(t100))
        # print(current)

        amplitudes, frequencies = extracter(x100, t100, num_peaks)
        # ⇓ replace with constitutive relationship between damping and current ⇓
        damping = 1.5*(current**2)
        # ⇑--------------------------------------------------------------------⇑
        power = 0
        for i in range(len(amplitudes)):
            power += (1/2)*damping*(amplitudes[i]**2)*((frequencies[i]*2*np.pi)**2)
        print(f"\nPower: {round(power, 4)}W")
        print(f"Peaks: amplitude of {np.round(amplitudes, 4)}m and frequency of {np.round(frequencies, 4)}Hz; Mean current: {current}A\n")
```

**Result of this code with two signals**

<img src="https://github.com/lgray14/Ocean-wave-energy-device/blob/main/images/power_ex.gif" height="500">


## GUI

After finishing the power extraction code, the next step was to integrate it with the existing graphical user interface (GUI) which plots the position of the buoy live. Initially, it plotted the live position and force reading from the force probe. My goal was to add the current to be plotted and the calculated power. Plotting the current turned out to be relatively trivial, since the current is read straight to the Arduino and can be passed directly to the GUI. I replaced the force reading with the current and just changed all of the labels to say current instead of force and all the appropriate units. At some point it would probably be useful to add the force probe functionality back to the Arduino readings and the GUI. From there the only major task was to add the power calculation functionality and have it plot. 

This proved to be fairly complex. I added in the FFT calculation functions, using Claude AI to help me sort through the GUI because I don't have a lot of experience reading functions and coding user interfaces -- Claude suggested using a buffer with the deque function from the collectons library in order to use a moving window for data collection. This makes it less computationally intensive than storing all of the position, time, and current values and only taking the last 100. With this and a few other minor modifications, I was able to implement live power readings into the GUI and print them to the bottom section where the raw position data is displayed.

Below is the code added for the FFT power analysis.

```python
def run_fft_power_analysis(self):
    if len(self.fft_buffer) < 20:
        return

    times, positions, currents = zip(*self.fft_buffer)
    # print(self.fft_buffer)
    times = np.array(times)
    positions = np.array(positions) / 100.0  # Convert cm to meters
    currents = np.array(currents) / 1000.0   # Convert mA to A
    # times -= times[0]  # Normalize time

   if len(times) >= 100:
        # Check time intervals
        time_window = times[-100:]
        time_intervals = np.diff(time_window)
        min_interval = np.min(time_intervals)
        max_interval = np.max(time_intervals)
        avg_interval = np.mean(time_intervals)
            
        if min_interval < 0.001:  # Less than 1ms between samples is suspicious
            self.report(f"Warning: Very small time interval detected: {min_interval:.6f}s")
            return
                
       if max_interval / min_interval > 10:  # Suspicious time interval variation
           self.report(f"Warning: Large time interval variation - min: {min_interval:.3f}s, max: {max_interval:.3f}s, avg: {avg_interval:.3f}s")
           return  # Skip this FFT calculation to avoid potential spikes
            
      # Create evenly spaced time points for interpolation
      desired_sample_rate = 50  # Hz
      desired_interval = 1.0 / desired_sample_rate
      t_uniform = np.linspace(time_window[0], time_window[-1], int((time_window[-1] - time_window[0]) / desired_interval))
            
      # Interpolate position data to uniform time points
      from scipy.interpolate import interp1d
      position_interp = interp1d(time_window, positions[-100:], kind='cubic', bounds_error=False)
      positions_uniform = position_interp(t_uniform)
            
      # For current, use linear interpolation since it might be noisier
      current_interp = interp1d(time_window, currents[-100:], kind='linear', bounds_error=False)
      currents_uniform = current_interp(t_uniform)
            
      # Use interpolated data for calculations
      mean_current = np.mean(currents_uniform)
      norm_times = t_uniform - t_uniform[0]
            
      # Debug interpolation quality
      if np.any(np.isnan(positions_uniform)):
          self.report("Warning: Interpolation produced NaN values")
          return
                
      # Apply Hanning window to reduce spectral leakage
      window = np.hanning(len(positions_uniform))
      positions_windowed = positions_uniform * window
            
      # Compensate for window amplitude reduction
      window_correction = .4/.294 # Hann window reduces amplitude -- experimentally determined correction factor
            
      amplitudes, frequencies = extracter(positions_windowed, norm_times, k=1)
      amplitudes = [amp * window_correction for amp in amplitudes]  # Apply correction

      if not amplitudes or not frequencies:
          self.report("FFT: No significant frequency detected.")
         return

      # Apply frequency stability check
      if hasattr(self, 'last_freq') and self.last_freq is not None:
          freq_change = abs(frequencies[0] - self.last_freq)
         if freq_change > 1:  # More than 1 Hz change
             self.report(f"Warning: Large frequency change: {freq_change:.2f} Hz")
             return
      self.last_freq = frequencies[0]

      damping = 1.5 * mean_current**2
      power = sum(0.5*damping*(amp**2)*(2*np.pi*f)**2 for amp, f in zip(amplitudes, frequencies))

      # Store the absolute time for the power value
      self.power_buffer.append((self.fft_buffer[-1][0], power))

      self.report(f"Power: {power:.4f} W | Freq: {frequencies[0]:.2f} Hz | Amp: {amplitudes[0]:.4f} m | Damping: {damping} N/m/s")
```

Adding the FFT-enabling functions to the GUI was not too hard, but adding the power to the current plot was a little tricky. I decided to put the plot the power on the same graph as the current, so I had to twin the axes so that there was a dependent variable axis for power output in Watts. In the end most of the edits had to happen in the *update_live()* function in the plotter library that enables the GUI. Once I was able to pass the power values to the *update_live()* function, I was able to plot them but I noticed a couple issues.

1. The power was spiking or reading crazy values when the position lagged slightly. It was tricking it into thinking the frequency went up really high or low when a lag meant that a bunch of position points were read all at once from the position sensor.
2. The plot was not showing a line, only individual points were appearing.
3. The plot was not shifting as time went on, instead the time axis was staying still. The GUI before automatically moved the time axis as new data came in.

To fix issue number 1, I added filtering to make sure that power was only being calculated and passed to *update_live()* when the time gaps were consistent -- not too large or too small. This happens in *run_fft_power_analysis()*, the function in the above code. This has done a good job of preventing major spikes and abnormalities, but it does shut off the power reading at any lag. I tried adding a timestamp to the position data as it comes in, but this caused massive lagging and a couple weird errors I couldn't sort out. This method makes the GUI run faster and more smoothly but in an ideal world the incoming data could be timestamped to avoid this issue entirely. To fix issue numbers 2 and 3, I was out of my depth and asked Claude for help as well. Together, we came up with the following block of code which fixed both problems. Using power_x and power_vals to plot and ensuring *self.line_power.set_visible* was set to *(True)* fixed the line issue, and using the aligned power times fixed the issue with the time axis not moving.
```python
if hasattr(self, 'line_power') and self.line_power is not None:
    # print("got here")
    if len(power_times) > 0:
        # Align power x-axis to match the current plot's time window
        if len(t_raw) > 0:
            t0 = t_raw[-1]
            power_x = -(power_times - t0)
            self.line_power.set_data(power_x, power_vals)
            # print("Power times (aligned)", power_x)
        else:
            self.line_power.set_data(power_times, power_vals)
        self.line_power.set_visible(True)
        # print('set data to t and power_vals')

        # Optionally, keep the power axis y-limits fixed or auto-adjust
        adjust_ylim(self.ax_current_twin, power_vals)
        # # Also, set the x-limits to match the current plot
        # self.ax_current_twin.set_xlim(self.time_window, 0.0)
    else:
        self.line_power.set_data((0, 0), (np.nan,) * 2)
        self.line_power.set_visible(False)
        # print('cleared data')
self.draw()
```
Once I had figured all that out, the GUI showed the position on the left plot, and the current and power on the right plot!

<img src="https://github.com/lgray14/Ocean-wave-energy-device/blob/main/images/sinusoid.png" height="500"> 

The very final step was to change the recording and saving data function to include the power. I did this by adding a column to the recording function, *toggle_recording*, for the power. Then, in the *save_data* function, I made sure the power is included and saved it with six decimal places. All of the other data save with two decimal places but since the power values are pretty low I wanted to make sure enough information about them were being saved. 


## Final

All of the relevant code I worked on this semester for the Arduino, python testing, and GUI are all housed in the 'code' folder of this GitHub. Wiring diagrams for all of the electronic components can be found online or on the linked informational webpages, and more information about the wave energy device and its design can be found in [Penelope Herrero-Marques' thesis](https://dspace.mit.edu/handle/1721.1/156654). 

The next steps for this project are to enable it to generate electricity. Having the ability to vary the damping is useful because in most real-life scenarios, the power takeoff (PTO) device generating electricity would have a fixed damping coefficient. Being able to vary the damping ratio in experiments means we can simulate a range of different PTO devices and determine what damping coefficient maximizes our system. Knowing this, we can proceed to create a device that matches this damping ratio. Another possibility would be to develop a device that can vary it's own damping ratio. This could be advantageous because in different scenarios, different coefficients could be more favorable, so having the ability to vary damping coefficient could make a generator more efficient. 

Two main options seem most reasonable for a PTO system. The first is a direct drive linear generator, where permanenent magnets would be installed on the moving element and pass through coils as the buoy bobbed. The relative motion would generate electricity in the coils. The cons of this approach are that there is inevitably some loss in the airgap, and the technology is not very advanced for a low velocity system like ours. It would likely be fabricated in house and be more expensive. However, there is no need for mechanical transformations with this approach, reducing mechanical losses in the electricity generated. The second option would be a combination of a linear to rotary conversion and a rotary generator. This would convert the linear motion of the bobbing device into a rotational motion, then use that rotational motion to drive a more traditional rotary generator. This option could be more cost effective and easier to build, but would have a lot of mechanical losses in the linear to rotary conversion. Since the system is small and produces a small amount of energy, this is less than ideal. Additionally, since the bobbing motion is pretty low velocity compared to, say, a piston in a traditional engine, it might be hard to find an effective way of converting the linear motion to rotational motion. Overall though, these would be two achievable options for achieving electricity generation with this system.

Additionally, characterizing the peak power conditions is another step to be taken before moving on to electricity generation. As of the end of the summer, testing had begun on the device to examine the power output at different wave frequencies and damping ratios. Analysis remains to be done on this data, spearheaded by Rafael Tejeda.
