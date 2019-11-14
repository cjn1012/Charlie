clear all; clc; close all; 

addpath('C:\Users\Student\OneDrive\College\Senior Year\S Lab\Lab 4\1');
addpath('C:\Users\Student\OneDrive\College\Senior Year\S Lab\Lab 4\2.1.3');
addpath('C:\Users\Student\OneDrive\College\Senior Year\S Lab\Lab 4\2.1.4');
addpath('C:\Users\Student\OneDrive\College\Senior Year\S Lab\Lab 4\2.2.1');
addpath('C:\Users\Student\OneDrive\College\Senior Year\S Lab\Lab 4\3.1');
addpath('C:\Users\Student\OneDrive\College\Senior Year\S Lab\Lab 4\3.2');

% addpath('C:\Users\User\Desktop\Lab 4\2.1.3');
% addpath('C:\Users\User\Desktop\Lab 4\1');
% addpath('C:\Users\User\Desktop\Lab 4\2.1.4');
% addpath('C:\Users\User\Desktop\Lab 4\3.1');
% addpath('C:\Users\User\Desktop\Lab 4\3.2');

%% Part 1: Potentiometer Accelerometer 
close all; 

acc_volt   = xlsread('Instrument_2 Capture 2019-10-31 15-14-13 Oscilloscope - Waveform Data.csv', 'A7:A100006');
t          = 0:1.024e-5:(1.024e-5) * 1e5;
t          = t(1:end-1);
dV         = 103/1000; 

% a.)
weight     = [10.1 60.1 110.1 160.1 210.1 260.1 310.1 360.1];
weight     = (weight./28.35);
increasing = [937 1015 1095 1173 1246 1330 1412 1484];
increasing = increasing./1000;
decreasing = [936 1015 1096 1175 1249 1331 1413 1491];
decreasing = decreasing./1000;

a          = polyfit(weight, increasing, 1);
sens       = a(1);
cal_curve  = a(1)*weight + a(2);

% b.) 
smooth               = 3;
[acc_volt, smooth]   = wsmooth(acc_volt, t, smooth);
acc_volt             = acc_volt - acc_volt(1);
min = -.1; 
for ii = 1:length(acc_volt)
    if acc_volt(ii) <= min
        min = acc_volt(ii);
    else
    end
end

percent_overshoot = 0.1524/0.6819;
zeta_Overshoot    = 0.45; 
meff              = (dV/a(1));

acc_volt          = acc_volt.*-1;

th = 0.01;
[peakLoc, peakMag] = peakfinder(acc_volt, th);
peakLoc(1)         = [];
peakMag(1)         = [];
% Damped Natural Frequency
n_peaks            = length(peakLoc);
Td                 = 1/(t(peakLoc(2)) - t(peakLoc(1)));
wd                 = Td*2*pi;                        % wd

figure()
plot(t, acc_volt, t(peakLoc), peakMag, 'd')
xlabel('Time [s]')
ylabel('Position')


for ii = 1:length(peakLoc)
    y(ii)          = log(peakMag(1)/peakMag(ii));
end

n = 0:length(y)-1;
damp_ratio         = zeros(1, length(peakLoc));

for jj = 1:length(peakLoc)
    num            = ((1/length(n))*log(peakMag(1)/peakMag(jj)));
    damp_ratio(jj) = num/(sqrt(4*pi^2 + num^2));
end

zeta               = mean(damp_ratio);
wn                 = wd/(sqrt(1 - zeta_Overshoot^2));    % wd = wn*sqrt(1-zeta^2) [rad/s]

spring_Const       = (meff*(wn/(2*pi))^2)/(32.2 * 12);   % ozf/in

% c.) 
sens1              = (meff*sens)/386;                    % V/(in/s2)

% d.) 
max_acc            = (increasing(end))/sens1;

% e.)
num                = [1/spring_Const];
den                = [meff ((2*zeta_Overshoot)/meff) spring_Const];
sys                = tf(num, den);
figure()
bode(sys)


% Plots 
figure()
plot(weight, increasing);
title('Output Voltage vs. Weight')
xlabel('Weight [oz_{f}]'); ylabel('e_{o} [V]')

figure() 
plot(t, acc_volt);
title('Output Voltage vs. Time')
xlabel('Time [s]'); ylabel('e_{o} [mV]')



%% Part 2 Piezoelectric Force Sensor 
% Part 2.1 

close all; 

piezo_w    = 0:200:1600;                        % grams
piezo_w    = (piezo_w./453.6);
voltage    = [0 23 44 69.8 85 106 126 141 167]; % mV 

b          = polyfit(piezo_w, voltage, 1);
sens2      = b(1); 
cal_curve2 = b(1)*piezo_w + b(2);

for ii = 1:length(voltage)
    
    residual(ii) = abs(voltage(ii) - cal_curve2(ii));
    
end

max_error  = max(residual);
percent_FS = max_error/voltage(end);

% Plots 
figure()
plot(piezo_w, voltage, 'k*-', piezo_w, cal_curve2, 'b*-') 
title('Calibration Curve for Piezoelectric Force Sensor');
xlabel('Weight [lb_{f}]'); ylabel('Voltage [mV]');

volt_data  = xlsread('Instrument Capture 2019-10-31 14-39-55 Oscilloscope - Waveform Data.csv', 'A6:A12008');
t          = 0:1e-6:(1e-6)*12002;
t          = t(1:end-1);

num_peaks  = 16; 
time_e     = 0.009177/.007161;

nat_freq   = num_peaks/time_e;     % Hz
nat_freq   = nat_freq*2*pi; 

figure()
plot(t, volt_data);
xlabel('Time [s]'); ylabel('Voltage [V]');
title('Voltage vs. Time for Piezoelectric Force Sensor');

% Part 2.2
impulse_mass  = 2; %lbf
impulse_volt  = xlsread('Impulse Loading.xlsx' ,'A7:A100007');
impulse_time  = t:2.4030e-5:(2.4030e-5)*100000;
impulse_time  = impulse_time(1:end-1);

plot(impulse_time, impulse_volt)

th = 0.005;
[peakLoc, peakMag] = peakfinder(impulse_volt, th);
peakLoc(1)         = [];
peakMag(1)         = [];
% Damped Natural Frequency
n_peaks            = length(peakLoc);
Td                 = 1/(impulse_time(peakLoc(2)) - impulse_time(peakLoc(1)));
wd                 = Td*2*pi;                        % wd
% 
figure()
plot(impulse_time, impulse_volt, impulse_time(peakLoc), peakMag, 'd')
xlabel('Time [s]')
ylabel('Voltage [mV]')


for ii = 1:length(peakLoc)
    y(ii)          = log(peakMag(1)/peakMag(ii));
end

n = 0:length(y)-1;
damp_ratio         = zeros(1, length(peakLoc));

for jj = 1:length(peakLoc)
    num            = ((1/length(n))*log(peakMag(1)/peakMag(jj)));
    damp_ratio(jj) = num/(sqrt(4*pi^2 + num^2));
end

impulse_zeta               = mean(damp_ratio);
impulse_wn                 = wd/(sqrt(1 - impulse_zeta ^2));  % 6/(1.6396 - 1.2787)  % this comes from counting 

impulse_springConst        = (impulse_mass*(impulse_wn/(2*pi))^2)/(32.2*12);
impulse_dampCoeff          = ((2*impulse_zeta)/impulse_wn); 

impulse_num                = [1/impulse_springConst];
impulse_den                = [impulse_mass/impulse_springConst impulse_dampCoeff/impulse_springConst 1];
impulse_sys                = tf(impulse_num, impulse_den);
[force, x]                 = impulse(impulse_sys, impulse_time(end));

figure() 
plot(x, (1/impulse_springConst)*force);%, impulse_time, impulse_volt)
xlabel('Time [s]'); ylabel('Voltage [mV]');
title('Simulated Impulse Response')





%% Part 3

close all; 
% a.)

lvt_volt  = xlsread('Instrument Capture 2019-10-31 14-50-31 Oscilloscope - Waveform Data.csv' ,'A6:A100008');
acce_volt = xlsread('Instrument Capture 2019-10-31 14-50-31 Oscilloscope - Waveform Data.csv' ,'B6:B100008');
time      = 0:1.168e-5:(1.168e-5)*100000;
time      = time(1:end-1);

figure()








