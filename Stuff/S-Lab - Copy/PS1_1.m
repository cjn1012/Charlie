clear all
close all
clc
%% Variables
C=1e-05; %Capacitance
R=1e04; %Resistance
e_i=5;  %Input voltage
num=100;
den=[1 30 100];
%% Setup Simulink
sim('PS1_sim.slx')
% figure% e_o=q2/C;
% plot(t,e_o)

% 
plot(t,e_0)
grid minor
title('Output Voltage Response to a Step Input, Q1 Part D')
xlabel('Time (sec)')
ylabel('Response (Volts)')