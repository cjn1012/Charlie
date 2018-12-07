clear all
close all
clc
%% Variables
Ri=10000;
Rf=1e6;
C=10e-6;

%% Transfer Function
num=(-Rf/Ri);
den=[Rf*C, 1];
trans=tf(num, den);
bode(trans)
grid minor
