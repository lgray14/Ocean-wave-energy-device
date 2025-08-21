// NEW SENSORS

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
#define calibration_factor 1250.0 //  -203700.0 //This value is obtained using the SparkFun_HX711_Calibration sketch
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
    // UNCOMMENT THE FOLLOWING IF GUI NOT WORKING - CURRENT WILL APPEAR INSTEAD OF FORCE
    Serial.print(0);
    Serial.print(' ');
    // 
    Serial.print(current_mA, 2); // current returns a float ==> GUI does not cope well, resolve
    Serial.print(' ');
    Serial.println(scale.get_units(), 2); //force scale.get_units() returns a float
  }
}
