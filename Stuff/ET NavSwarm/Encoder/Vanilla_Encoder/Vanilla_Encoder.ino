/*     Arduino Rotary Encoder Tutorial
 *      
 *  by Dejan Nedelkovski, www.HowToMechatronics.com
 *  
 */
 
 #define outputA 11
 #define outputB 7

 int counter = 0; 
 int aState;
 int aLastState;  
 int angle = 0;
 void setup() { 
   pinMode (outputA,INPUT);
   pinMode (outputB,INPUT);
   
   Serial.begin (9600);
   // Reads the initial state of the outputA
   aLastState = digitalRead(outputA);   
 } 

 void loop() { 
   aState = digitalRead(outputA); // Reads the "current" state of the outputA
   // If the previous and the current state of the outputA are different, that means a Pulse has occured
   if (aState != aLastState){     
     // If the outputB state is different to the outputA state, that means the encoder is rotating clockwise
     if (digitalRead(outputB) != aState) { 
       counter ++;
     } 
     
     //else {
      // counter --;
     //}
    
     Serial.println(" ");
     Serial.print("Last A: ");
     Serial.println(aLastState);
     Serial.print("Position: ");
     Serial.println(counter);
     Serial.print("Angle: ");
     angle = (counter*5.625);
     if (angle >= 360){
      angle = 0;
     }
     Serial.println(angle);
     Serial.print("A: ");
     Serial.println(digitalRead(outputA));
     Serial.print("B: ");
     Serial.println(digitalRead(outputB));
   } 
   aLastState = aState; // Updates the previous state of the outputA with the current state
 }
