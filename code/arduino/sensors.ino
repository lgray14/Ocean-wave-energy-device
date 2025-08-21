/*
  ReadAnalogVoltage

  Reads an analog input on pin 0, converts it to voltage, and prints the result to the Serial Monitor.
  Graphical representation is available using Serial Plotter (Tools > Serial Plotter menu).
  Attach the center pin of a potentiometer to pin A0, and the outside pins to +5V and ground.

  This example code is in the public domain.

  https://www.arduino.cc/en/Tutorial/BuiltInExamples/ReadAnalogVoltage
*/

#include "HX711.h"

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
}

// the loop routine runs over and over again forever:
void loop() {
  // read the input on analog pin 0:
  int sensorValue = analogRead(A0);
  int sensorValue2 = analogRead(A2);
  // Convert the analog reading (which goes from 0 - 1023) to a voltage (0 - 5V):
  // float voltage = sensorValue * (5.0 / 1023.0);
  // print out the value you read:
  // Serial.println(voltage);

  Serial.print(sensorValue);
  Serial.print(' ');
  Serial.print(sensorValue2);
  Serial.print(' ');
  Serial.println(scale.get_units(), 2); //scale.get_units() returns a float

  // Serial.println(sensorValue);

  // char buffer[40];

  // sprintf(buffer, "%4d %4d asoho", sensorValue, sensorValue2);
  // sprintf(buffer, "%4d %4d  23.1", sensorValue, sensorValue2);
  // Serial.println(buffer);
}
