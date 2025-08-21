void setup() {
  // put your setup code here, to run once:
  Serial.begin(115200);
}

float t = 0;
float amp = 350; // mm to simulate actual data read
float freq = .4;

float amp2 = 50;
float freq2 = 3;
void loop() {
  // put your main code here, to run repeatedly
  while(t <= 5){
    float dist = amp*sin(2*PI*freq*t); //+ amp2*sin(2*PI*freq2*t);
    // Serial.println(amp2*sin(2*PI*freq2*t));
    int rDist = round(dist);
    t += .1;
    Serial.print(rDist);
    Serial.print(' ');
    int secondDist = 0; //random(-50, 50);
    Serial.print(secondDist);
    Serial.print(' ');
    float current = 0;
    Serial.println(current);
    delay(100);
  }
  while(t > 5){
    float dist = amp*sin(2*PI*freq*t); // + amp2*sin(2*PI*freq2*t);
    int rDist = round(dist);
    t += .1;
    Serial.print(rDist);
    Serial.print(' ');
    int secondDist = 0; //random(-50, 50);
    Serial.print(secondDist);
    Serial.print(' ');
    float current = random(1450, 1550);
    Serial.println(current);
    delay(100);
  }
}
