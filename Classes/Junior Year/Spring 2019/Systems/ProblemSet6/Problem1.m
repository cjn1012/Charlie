%% Problem Set 5
clear all
close all
%% Problem 1 MATLAB
% System Constants


Numer = [1 1000 0];
Denom = [1 100010 (10^10+10^6) 10^11];

TF1 = tf(Numer,Denom);
figure(1)
bode(TF1)

Numer2 = [1 1001 10^4+10^3 10^4]
Denom2 = [1 10^8 10^6 0]

TF2 = tf(Numer2,Denom2);
figure(2)
bode(TF2)


