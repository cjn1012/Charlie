clear all
close all
clc
%% Declare Variables
E1=108e9;
E2=24.2e9;
E3=18.6e9;
v12=0.32;
v21=(v12*E2)/E1;
v23=0.18;
v32=(E3*v23)/E2;
v31=0.18;
v13=(E1*v31)/E3;
sig1=0;
T6=1200;
T5=900;
sig2=2018;
T4=-1400;
sig3=1600;
G23=9.2e9;
G31=12.4e9;
G12=16.5e9;
A= [1/E1 -v21/E2 -v31/E3 0 0 0
    -v12/E1 1/E2 -v32/E3 0 0 0
    -v13/E1 -v23/E2 1/E3 0 0 0
    0 0 0 1/G23 0 0
    0 0 0 0 1/G31 0
    0 0 0 0 0 1/G12];
b= [sig1
    sig2
    sig3
    T4
    T5
    T6];
b=b*6894.76; %Convert to Pa
strain=A*b

