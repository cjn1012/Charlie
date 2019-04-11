
clear all
close all
clc

ResX = 100;
ResY = 100;

T_0 = 20; % C
T_1 = 100; % C
L1 = 1; % Centimeters
L2 = .70; % cmeters
D = .50; % height cm
k1 = 2; % W/cmK
k2 = .50;  % W/cmK

%% call first function

[A, B, cc, rr] = func1(L1,L2, D, T_0, T_1,k1,k2,ResX,ResY);

%% call second function
[x] = func2(A, B);

bb = reshape(x, cc, rr)';
% h = A\NodeType

imagesc(bb)
colorbar
xlabel('L (cm)')
ylabel('D (cm)')
