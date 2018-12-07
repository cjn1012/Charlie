#define outputA 7
#define outputB 11

 int counter = 0; 
 int aState;
 int bState;
 int aLastState;  
 int angle = 0;
 void setup() { 
   pinMode (outputA,INPUT);
   pinMode (outputB,INPUT);
   
   Serial.begin (38400);
   // Reads the initial state of the outputA
   aLastState = digitalRead(outputA);   
 } 

 void loop() { 
   aState = digitalRead(outputA); // Reads the "current" state of the outputA
   // If the previous and the current state of the outputA are different, that means a Pulse has occured
      bState = digitalRead(outputB);
   if (aState != aLastState){     
     // If the outputB state is different to the outputA state, that means the encoder is rotating clockwise
     if (digitalRead(outputB) != aState) { 
       counter ++;
     } else {
       counter --;
     }
     Serial.print("Count: ");
     Serial.println(counter);
     Serial.print("Angle: ");
     angle = (1/64)*counter;
     if (angle >= 360){
      angle = 0;
     }
     Serial.println(angle);
   } 

 /*    Serial.print("B: ");
     Serial.println(bState);
     Serial.print("A: ");
     Serial.println(aState);
     */
   aLastState = aState; // Updates the previous state of the outputA with the current state
   //bLastState = bState;
 }
