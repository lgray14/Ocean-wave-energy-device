# Ocean wave energy device

* [Background](#Background)
* [Power extraction](#Power)
* [GUI and interface](#GUI)

### Background

This project is to design and test a wave energy converter that will use the motion of waves in the ocean to produce electrical energy. Ocean waves are a huge and largely untapped renewable energy source – wave energy could provide an estimated easily exploitable 500 gigawatts according to [CorPower](http://corpowerocean.com/a-short-history-of-wave-energy/) – and a major advantage of wave energy is that it is not dependent on weather conditions or time of day. Energy sources such as solar power and wind energy rely on variable conditions; solar depends on availability of sunlight, and wind depends on whether and how fast the wind is blowing. In contrast, in many coastal environments waves are constant throughout the day and night, meaning that power can be constantly generated and does not have to be stored to be distributed as much as other renewable energy forms. Wave energy thus has the capacity to provide a stable, reliable source of energy at scale. 

Building on the work done by The Vortical Flow Research Laboratory, this project tests and characterizes an existing wave-energy converter device. The type of wave energy device being studied is a point absorber, which consists of a floating buoy and a fixed rig. The relative motion of the two would drive a power takeoff device (PTO). In our system, the PTO is simulated by a linear damper - the damper dissipates the energy of the motion of the buoy using electromagnetic resistance. A standard PTO would have a fixed damping coefficient – this setup allows us to test multiple different damping coefficients, simulating different PTO options, without having to reconfigure the device. The goal is to determine the optimal configuration that produces the most power for this particular setup.

<img src="https://github.com/lgray14/Ocean-wave-energy-device/blob/main/images/buoy.jpg" height="300">

<sup> *The device -- a modified Home Depot bucket which has been constrained to move vertically* </sup>


My role in this project has been to run tests, conduct analysis of the data collected, and to add functionality to the device. By adding sensors and writing code, I have been able to calculate and display the instantaneous power output that the device is producing. Through lab tests, I have helped to characterize the linear damper to get an exact relationship between the current to the electromagnets and the damping ratio, and have verified that the code I have worked on functions in lab scenarios. Furthermore, I have integrated my work with the Graphical User Interface (GUI) to better communicate the results. 


### Power


### GUI
