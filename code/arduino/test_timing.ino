// SAMPLING RATE TEST

#include "Simpletimer.h"

Simpletimer timer1{};

void setup() {
  // put your setup code here, to run once:
  Serial.begin(115200);
}


float sampling_rate = 20; //Hz
float gap = (1/sampling_rate)*1000;

void loop() {
  if (timer1.timer(gap)) {
    Serial.print("Current time: ");
    Serial.println(millis());
  }
}