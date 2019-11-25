clear all
close all

%% Reading in Data for c,d,e
P_TachOut = xlsread('P-Control.csv', 'A7:A100007');
P_MotorIn = xlsread('P-Control.csv', 'B7:B100007');
P_Time    = xlsread('P-Control.csv', 'C7:C100007');
figure(1)
plot(P_Time,P_TachOut)
hold on
plot(P_Time,P_MotorIn)
legend('Tach Out','Motor In')

I_TachOut = xlsread('I-Control.csv', 'A7:A100007');
I_MotorIn = xlsread('I-Control.csv', 'B7:B100007');
I_Time    = xlsread('I-Control.csv', 'C7:C100007');
figure(2)
plot(I_Time,I_TachOut)
hold on
plot(I_Time,I_MotorIn)
legend('Tach Out','Motor In')

PI_TachOut = xlsread('PI-Control.csv', 'A7:A100007');
PI_MotorIn = xlsread('PI-Control.csv', 'B7:B100007');
PI_Time    = xlsread('PI-Control.csv', 'C7:C100007');
figure(3)
plot(PI_Time,PI_TachOut)
hold on
plot(PI_Time,PI_MotorIn)
legend('Tach Out','Motor In')

%% a

% drawing the block diagrams of each controller attached in  folder


%% b 
% Root Locus Plots for each controller

Ri  = 100000;  %ohm
Rf  = 150000;  %ohm
C1  = .033e-6; %farad
C2  = .010e-6; %farad
K   = 2;   %unitless
tau = 0.0103;  %sec
K_m = 1; %changing gain []

P_sys  = tf([K_m*Rf*K],[tau 1]);
I_sys  = tf([K_m*K],[tau*C1 C1 0]);
PI_sys = tf([K_m*K*C2*Rf K_m*K],[tau*C2 C2 0]);

figure(4)
rlocus(P_sys)
figure(5)
rlocus(I_sys)
figure(6)
rlocus(PI_sys)

%% c - finding roots

[rp,kp] = rlocus(P_sys);
Proot1_p = rp(1);
Proot1_k = kp(1);

[rp,kp] = rlocus(I_sys);
Iroot1_p = rp(1,1);
Iroot1_k = kp(1);
Iroot2_p = rp(2,1);
Iroot2_k = kp(1);

[rp,kp] = rlocus(PI_sys);
PIroot1_p = rp(1,1);
PIroot1_k = kp(1);
PIroot2_p = rp(2,1);
PIroot2_k = kp(1);






