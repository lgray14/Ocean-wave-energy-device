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

Building on the work done by The Vortical Flow Research Laboratory, this project tests and characterizes an existing wave-energy converter device. The type of wave energy device being studied is a point absorber, which consists of a floating buoy and a fixed rig. The relative motion of the two would drive a power takeoff device (PTO). In our system, the PTO is simulated by a linear damper - the damper dissipates the energy of the motion of the buoy using electromagnetic resistance. A standard PTO would have a fixed damping coefficient – this setup allows us to test multiple different damping coefficients, simulating different PTO options, without having to reconfigure the device. The goal is to determine the optimal configuration that produces the most power for this particular setup.

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

<img src="https://github.com/lgray14/Ocean-wave-energy-device/blob/main/images/buoy.jpg" height="300">

<sup> *When the bucket moves in the tank, its position varies sinusoidally as shown here* </sup>


### Sensors
To get the position data, a position sensor must be equipped with the device. The device previously had a distance sensor that worked by turning as the buoy bobbed up and down -- the contact with the frame made it turn and the number of turns corresponded to the distance travelled. However, this sensor added a lot of friction to the system and we decided it was better to use a [time of flight (ToF) sensor](https://www.adafruit.com/product/3317). The VL53L0X ToF sensor seemed like a good choice, since the range need was fairly short and the accuracy need pretty high, and the sensor itself is small and lightweight. The sensor can take readings at different frequencies, so it needed to be tuned to take as many readings as possible in a certain band of accuracy. After some trial, the following code represents the best balance between speed of readings and accuracy:
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

I placed the sensor on the frame of the wave device to be able to measure the distance to the moving buoy. Having the sensor fixed instead of moving with the system made the wiring less complicated and more secure.

<img src="https://github.com/lgray14/Ocean-wave-energy-device/blob/main/images/buoy.jpg" height="300">


### FFT

### Filtering


## GUI
