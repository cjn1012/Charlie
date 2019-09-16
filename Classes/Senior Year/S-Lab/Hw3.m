clear all
close all

% Prob 1b

R = 1000
C = .00001

figure(1)
H = tf([-R*C,0],[1])
bode(H)

[Z,gain] = zero(H)

