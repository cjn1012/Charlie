% Charlie Nitschelm

% clear all
% close all
clc
%% This Script Pulls from the MatrixMaker and MatrixSolver Functions and displays a steady 2D temperature figure

% Input Parameters
ResX = 80;
ResY = 30;
T_0 = 100+273.15; % C
T_1 = 200+273.15; % C
L1 = 1; % Meters
L2 = .70; % Meters
D = .50; % Meters
k1 = 200; % W/mK
k2 = 100;  % W/mK

%% call first function

[A, B, cc, rr,NodeType] = MatrixMaker(L1,L2, D, T_0, T_1,k1,k2,ResX,ResY);

%% call second function
[x] = AB(A, B);

bb = reshape(x, cc, rr)';

figure(1)
imagesc([0,L1+L2],[0,D],bb)
colorbar
xlabel('Length (m)')
ylabel('Depth (m)')
