clear all
close all

num=1;
den=[1/(2000^2), 2.0000e-05, 1];
tf=tf(num,den);
bode(tf)
title('Problem 4')
grid minor

text(3*10^3, 150, '\omega_n = 2*10^2 rad/s')
text(2*10^2, -45, 'dB/dec after \omega_n = -40 dB/dec')