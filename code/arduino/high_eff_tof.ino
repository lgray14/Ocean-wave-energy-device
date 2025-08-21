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
  // if (Status == VL53L0_ERROR_NONE) {
  //   Status = VL53L0_SetLimitCheckValue(pMyDevice,  VL53L0_CHECKENABLE_SIGNAL_RATE_FINAL_RANGE, (FixPoint1616_t)(0.25*65536));
  // }
  // if (Status == VL53L0_ERROR_NONE) {
  //   Status = VL53L0_SetLimitCheckValue(pMyDevice, VL53L0_CHECKENABLE_SIGMA_FINAL_RANGE, (FixPoint1616_t)(18*65536));
  // }
  // if (Status == VL53L0_ERROR_NONE) {
  //   Status = VL53L0_SetMeasurementTimingBudgetMicroSeconds(pMyDevice, 200000);
  // }

  lox.setMeasurementTimingBudgetMicroSeconds(100000);
}

void loop() {
  VL53L0X_RangingMeasurementData_t measure;
  lox.rangingTest(&measure, false); // pass in 'true' to get debug data printout!

  int distMeasured = 0;
  if (measure.RangeStatus != 4) {  // phase failures have incorrect data
    distMeasured = measure.RangeMilliMeter;
  }  

  // Serial.print('Distance: ');
  Serial.println(distMeasured);
  // Serial.print(' mm    Time:');
  // Serial.println(millis());
  // Serial.println(' s');
  Serial.println(' ');
}
