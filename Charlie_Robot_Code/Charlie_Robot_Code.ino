// Sweep


#include <Servo.h>

Servo myservo;  // create servo object to control a servo
Servo myservo2; // twelve servo objects can be created on most boards

int pos  = 0;    // variable to store the servo position
int pos2 = 0;


void setup() {
myservo.attach(9);    // attaches the servo on pin 9 to the servo object elbow
myservo2.attach(10);  // attached the servo on pin 10 to the servo object hand
/*up is 50 down is 150
 * open is 50
 * clamp is 90
 */
}

void loop() {
// Up and Open

for (pos=150; pos>50;pos-=1) {
  myservo.write(pos);
  delay(10) ; 
}

delay(1000);
myservo2.write(50);
delay(1000);

for (pos=50; pos<150; pos+=1) {
  myservo.write(pos);
  delay(10);
}
delay(1000);

// Closed
myservo2.write(90);
delay(1000);

    

}
