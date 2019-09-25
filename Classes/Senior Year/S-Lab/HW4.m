clear all
close all

ks = 4e-7;
t1 = 1/6000;
t2 = 1/60000;

sys = tf([ks,0],[t1*t2,t1+t2,1])

bode(sys)