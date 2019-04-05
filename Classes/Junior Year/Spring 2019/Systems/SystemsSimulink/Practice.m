% Problem Set 5 --- Problem 5

clear all;
close all;


% Mass Spring Damper

M = 1;
B = 1;
k = 10;

% Plotting Step Response of X(s) and sX(s)

% Transfer Function

% Define Numerator
num = [1];

% Define Denominator
den = [M,B,k]; % Represents polynomial, want to make transferfunction!

G_1 = tf(num,den)

figure(1)
step(G_1) % Creates step response for transfer function, assumes unit step input


% Next G_2

num2 = [1,0];

G_2 = tf(num2,den)
figure(2)
stepinfo(G_2)
grid



A = [0,1;-k/M,-B/M];

B = [0;1/M];

C = [10;.01];

D = [0;0];

sys = ss(A,B,C,D)

ss2tf
tf2ss











