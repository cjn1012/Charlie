clear all; close all;
% Filters
Voltage = [0,-.1,-.2,-.3,-.4,-.5,-.6,-.7,-.8,-.9,-1,-1.1,-1.2,-1.3,-1.4,-1.5,-1.6,-1.7,-1.8,-1.9,-2];
DarkCurrent = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
Current577 = [49,33,20,10,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0].*10^-11;
Current435 = [35,33,31.5,29.5,27.5,25.5,23.5,21,19.5,17,15,13.5,12,10.5,9,7,5,3.3,2.1,1,0]*10^-10;
Current365 = [82,72,64,58,50,43,38,30,24,18,13,8,4,.6,.33,.13,.06,0,0,0,0]*10^-10;
Current546 = [26,20,12,5,.84,.08,.012,.008,.008,.008,.008,.008,.008,.008,.008,.008,.008,.008,.008,.008,.008]*10^-10;
Current404 = [73,63,53,43,37,29,21,15,9,4,.51,.22,.11,-.7,-.12,-.14,-.15,-.16,-.16,-.16,-.16]*10^-10;

figure(1)
plot(Voltage,Current365)
hold on
plot(Voltage,DarkCurrent)
legend('365nm','Dark Current')

figure(2)
plot(Voltage,Current404)
hold on
plot(Voltage,DarkCurrent)
legend('404nm','Dark Current')

figure(3)
plot(Voltage,Current435)
hold on
plot(Voltage,DarkCurrent)
legend('435nm','Dark Current')

figure(4)
plot(Voltage,Current546)
hold on
plot(Voltage,DarkCurrent)
legend('546nm','Dark Current')

figure(5)
plot(Voltage,Current577)
hold on
plot(Voltage,DarkCurrent)
legend('577nm','Dark Current')

figure(6)
plot(Voltage,Current365)
hold on
plot(Voltage,Current404)
plot(Voltage,Current435)
plot(Voltage,Current546)
plot(Voltage,Current577)
legend('365nm','404nm','435nm','546nm','577nm','location','northwest')
















