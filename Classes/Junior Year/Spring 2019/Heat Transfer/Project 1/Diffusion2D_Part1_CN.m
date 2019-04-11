clear all
close all
clc
%% This Script Pulls from the MatrixMaker and MatrixSolver Functions and displays a steady 2D temperature figure

% Input Parameters
ResX = 30;
ResY = 30;
T_0 = 20; % C
T_1 = 100; % C
L1 = 1; % Meters
L2 = .70; % Meters
D = .50; % Meters
k1 = 200; % W/mK
k2 = 50;  % W/mK

%% call first function

[A, B, cc, rr] = func1(L1,L2, D, T_0, T_1,k1,k2,ResX,ResY);

%% call second function
[x] = func2(A, B);

bb = reshape(x, cc, rr)';

figure(1)
imagesc([0,L1+L2],[0,D],bb)
colorbar
xlabel('Length (m)')
ylabel('Depth (m)')
