/*
  ReadAnalogVoltage

  Reads an analog input on pin 0, converts it to voltage, and prints the result to the Serial Monitor.
  Graphical representation is available using Serial Plotter (Tools > Serial Plotter menu).
  Attach the center pin of a potentiometer to pin A0, and the outside pins to +5V and ground.

  This example code is in the public domain.

  https://www.arduino.cc/en/Tutorial/BuiltInExamples/ReadAnalogVoltage
*/

#include "HX711.h"
#include <Wire.h>
#include <Adafruit_INA219.h>
#include <Adafruit_VL53L0X.h>

//hall effect sensor
Adafruit_INA219 ina219;
// ToF sensor
Adafruit_VL53L0X lox = Adafruit_VL53L0X();

#define calibration_factor 1250.0 //  -203700.0 //This value is obtained using the SparkFun_HX711_Calibration sketch

#define DOUT  3
#define CLK  2

HX711 scale;


// the setup routine runs once when you press reset:
void setup() {
  // initialize serial communication at 9600 bits per second:
  Serial.begin(115200);

  // Serial.available();

  scale.begin(DOUT, CLK);
  scale.set_scale(calibration_factor); //This value is obtained by using the SparkFun_HX711_Calibration sketch
  scale.tare(); //Assuming there is no weight on the scale at start up, reset the scale to 0

  Serial.println("Readings:");

  // initialize Hall Effect sensor
  ina219.begin();
  // set up ToF sensor
  if (!lox.begin()) {
    Serial.println(F("Failed to boot VL53L0X"));
    while(1);
  }
}

// the loop routine runs over and over again forever:
void loop() {
  // read the input on analog pin 0:
  // int sensorValue = analogRead(A0);

// CANNOT TELL WHAT THIS ONE DOES //
  // int sensorValue2 = analogRead(A2);

  // Convert the analog reading (which goes from 0 - 1023) to a voltage (0 - 5V):
  // float voltage = sensorValue * (5.0 / 1023.0);
  // print out the value you read:
  // Serial.println(voltage);

  // new ToF sensor
  VL53L0X_RangingMeasurementData_t measure;

  lox.rangingTest(&measure, false); // pass in 'true' to get debug data printout!

  
    int distMeasured = 0;
    if (measure.RangeStatus != 4) {  // phase failures have incorrect data
      distMeasured = measure.RangeMilliMeter;
    }

    // Hall effect sensor
    float current_mA = 0;
    current_mA = ina219.getCurrent_mA();

    Serial.print(distMeasured);
    Serial.print(' ');
    // Serial.print(sensorValue2);
    // Serial.print(' ');
    Serial.print(0); // print hall effect sensor results
    Serial.print(' ');
    // ALSO DK WHAT THIS DOES //
    Serial.println(current_mA, 2); //scale.get_units() returns a float
  // }
  // else {
  //   Serial.println(" out of range "); // if too far away, print 'out of range'
  // }
    
  // delay(10);

  // Serial.println(sensorValue);

  // char buffer[40];

  // sprintf(buffer, "%4d %4d asoho", sensorValue, sensorValue2);
  // sprintf(buffer, "%4d %4d  23.1", sensorValue, sensorValue2);
  // Serial.println(buffer);
}
