%% Problem Set 5
clear all
close all
%% Problem 1 MATLAB
% System Constants
k = 3;
J = 3;
B = [3,6,12,0]; % Under, critically,over,undamped
K = 1;

% Calcs
Omega_n = sqrt(k/J);
Zeta = B./(2*sqrt(k*J));

Numer = K;
Denom = [1./Omega_n^2,2.*Zeta(1)/Omega_n,1;...
    1./Omega_n^2,2.*Zeta(2)/Omega_n,1;...
    1./Omega_n^2,2.*Zeta(3)/Omega_n,1;...
    1./Omega_n^2,2.*Zeta(4)/Omega_n,1];

TF1 = tf(Numer,Denom(1,:));
SS1 = ss(TF1);
TF2 = tf(Numer,Denom(2,:));
SS2 = ss(TF2);
TF3 = tf(Numer,Denom(3,:));
SS3 = ss(TF3);
TF4 = tf(Numer,Denom(4,:));
SS4 = ss(TF4);

% Vectors
t = linspace(0,20,2000);
u = ones(1,length(t));


figure(1)
lsim(SS1,u,linspace(0,16,2000),[1 1])
grid
title('Underdamped System - J=3,k=3,b=3')
ylabel('Amplitude')
xlabel('Time (s)')

figure(2)
lsim(SS2,u,linspace(0,12,2000),[1 1])
grid
title('Critically Damped System - J=3,k=3,b=6')
ylabel('Amplitude')
xlabel('Time (s)')

figure(3)
lsim(SS3,u,linspace(0,20,2000),[1 1])
grid
title('Overdamped System - J=3,k=3,b=12')
ylabel('Amplitude')
xlabel('Time (s)')

figure(4)
lsim(SS4,u,linspace(0,8,2000),[1 1])
grid
title('Undamped System - J=3,k=3,b=0')
ylabel('Amplitude')
xlabel('Time (s)')





