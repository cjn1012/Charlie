clear all; close all; clc; 

% addpath('C:\Users\User\Desktop\Charlie\Classes\Senior Year\S-Lab\Lab3\1.2');
% addpath('C:\Users\User\Desktop\Charlie\Classes\Senior Year\S-Lab\Lab3\2.1.2');
% addpath('C:\Users\User\Desktop\Charlie\Classes\Senior Year\S-Lab\Lab3\2.1.4');
% addpath('C:\Users\User\Desktop\Charlie\Classes\Senior Year\S-Lab\Lab3\2.2.4');

% On SEDS Room Comp
addpath('C:\Users\User\Desktop\Charlie\Classes\Senior Year\S-Lab\Lab3\1.2');
addpath('C:\Users\User\Desktop\Charlie\Classes\Senior Year\S-Lab\Lab3\2.1.2');
addpath('C:\Users\User\Desktop\Charlie\Classes\Senior Year\S-Lab\Lab3\2.1.4');
addpath('C:\Users\User\Desktop\Charlie\Classes\Senior Year\S-Lab\Lab3\2.2.4');


%% Linear Velocity Transducer (LVT)
close all; 

t   = 0:1.201e-5:(1.201e-5)*100000;
t   = t(1:end-1);

LVT = xlsread('Instrument Capture 2019-10-16 14-22-52 Oscilloscope - Waveform Data.csv', 'A6:A100006');

figure()
hold on
plot(t, LVT);
xlabel('Time [s]'); ylabel('LVT Voltage [V]');
xlim([t(1) t(end)]); grid on;
title('LVT Voltage vs. Time');

%e.) 
x_rest      = 2.29;          %inches
x_raised    = 0.70;          %inches
x_dist      = x_rest - x_raised; 

ii = 1; 
while LVT(ii) < 0.03 && LVT(ii) > -0.05
    ii = ii + 1;
end
start = ii;

min = -3;
for ii = 1:length(LVT)
    if(LVT(ii) < min)
        min = LVT(ii);
        endevent = ii;
    end
end

factor = 4000;
plot(t(start), LVT(start), '*', t(start+factor), LVT(start+factor), '*');
hold off

lvt_sens = ((LVT(start+factor) - LVT(start))/(t(start+factor) - t(start)))/386; %V/(in/s)

lvt_vel = LVT./lvt_sens;

lvt_pos = cumtrapz(t, lvt_vel);


figure()
hold on 
title('Core Velocity vs. Time')
xlabel('Time [s]'); ylabel('Core Velocity [in/s]');
xlim([t(1) t(end)])
plot(t, lvt_vel); grid on; hold off;

figure()
plot(t, lvt_pos)
title('Core Position vs. Time');
xlabel('Time [s]'); ylabel('Core Position [inch]');
xlim([t(1) t(end)])
grid on; 

% g.) Using the log decrement process I found the damped natural frequency
% and natural frequency. For the spring constant I took advantage of the
% fact that the step response for position is second order so the
% steady-state value is equal to 1/K (this is a mass-spring-damper
% system.) Since they wanted it in lbf/in, I used an online converter from
% grams to lbm and then multiplied by the gravitational constant.
th = 0.01;
[peakLoc, peakMag] = peakfinder((lvt_pos - mean(lvt_pos(end-1000:end))), th);
peakLoc(1)         = [];
peakMag(1)         = [];

% Damped Natural Frequency
n_peaks            = length(peakLoc);
Td                 = 1/(t(peakLoc(2)) - t(peakLoc(1)));
wd                 = Td*2*pi;                        % wd

% figure()
% plot(t, (lvt_pos - mean(lvt_pos(end-1000:end))), t(peakLoc), peakMag, 'd')
% xlabel('Time [s]')
% ylabel('Position')


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
wn                 = wd/(sqrt(1 - zeta^2));    % wd = wn*sqrt(1-zeta^2)
amountcompressed = 1.59-lvt_pos(end);
lvt_K              = (1/amountcompressed)*(75.6/453.6); 

% The spring constant is 3.5 lbf/in. This seems to make sense. The material
% is foam, so for every 2.5 pounds applied the foam depresses 1 inch. The wn
% also seems reasonable. If you count the peaks you get just about 4
% peaks/0.5 seconds. If you divide wn by 2pi you get about 7 Hz or 7
% oscillations a second. These values are close. 


% Here i am finding the point when mass hits the foam, and then I will multiply by b and mass and thats the initial force

ss = lvt_pos(end);
%x = (75.6/453.592)*32.2/(lvt_K*12)
foam_begin = ss-amountcompressed;

for i = 1:length(lvt_pos)
    if lvt_pos(i) > foam_begin
        index = i;
        break
    end
end

vel_foam = lvt_vel(index);

force = vel_foam*2*zeta/wn


%% LVDT 
close all; clear all;

% addpath('C:\Users\Lucas\OneDrive\College\Senior Year\S Lab\Lab 3\2.1.2\above null');
% addpath('C:\Users\Lucas\OneDrive\College\Senior Year\S Lab\Lab 3\2.1.2\null');
% addpath('C:\Users\Lucas\OneDrive\College\Senior Year\S Lab\Lab 3\2.1.4\above null');
% addpath('C:\Users\Lucas\OneDrive\College\Senior Year\S Lab\Lab 3\2.2.4');

% SEDS Room Computer 
addpath('C:\Users\User\Desktop\Charlie\Classes\Senior Year\S-Lab\Lab3\2.1.2\above null');
addpath('C:\Users\User\Desktop\Charlie\Classes\Senior Year\S-Lab\Lab3\2.1.2\null');
addpath('C:\Users\User\Desktop\Charlie\Classes\Senior Year\S-Lab\Lab3\2.1.4\above null');
addpath('C:\Users\User\Desktop\Charlie\Classes\Senior Year\S-Lab\Lab3\2.2.4');

rawNull_coil1 = xlsread('Instrument Capture 2019-10-16 14-43-38 Oscilloscope - Waveform Data.csv', 'A6:A6851');
rawNull_coil2 = xlsread('Instrument Capture 2019-10-16 14-43-38 Oscilloscope - Waveform Data.csv', 'B6:B6851');
raw_coil1 = xlsread('Instrument Capture 2019-10-16 14-44-40 Oscilloscope - Waveform Data.csv', 'A6:A6851');
raw_coil2 = xlsread('Instrument Capture 2019-10-16 14-44-40 Oscilloscope - Waveform Data.csv', 'B6:B6851');

filter_coil = xlsread('Instrument Capture 2019-10-16 14-52-59 Oscilloscope - Waveform Data.csv', 'A6:A6851');
mod_coil = xlsread('Instrument Capture 2019-10-16 14-52-59 Oscilloscope - Waveform Data.csv', 'B6:B6851');

coil_response = xlsread('Instrument Capture 2019-10-16 15-26-56 Oscilloscope - Waveform Data.csv', 'A6:A100006'); 


t = 0:3.8e-7:(6844*(3.8e-7));
t = t(1:end-1);

tResponse = 0:5.721e-5:(1e5)*5.721e-5;
tResponse = tResponse(1:end-1);


figure()
subplot(2,1,2)
plot(t, raw_coil1, t, raw_coil2);
title('Above Null Raw Coil Outputs of LVDT');
ylabel('Voltage [V]'); xlabel('Time [s]');
xlim([t(1) t(end)]); 
legend('Coil 1', 'Coil 2'); grid on;

subplot(2,1,1)
plot(t, rawNull_coil1, t, rawNull_coil2);
title('Null Coil Raw Outputs of LVDT');
ylabel('Voltage [V]'); 
xlim([t(1) t(end)]);
legend('Coil 1', 'Coil 2'); grid on;

figure()
plot(t, mod_coil);
title('Coil Output for LVDT After Demodulation');
xlabel('Time [s]'); ylabel('Voltage [V]');
xlim([t(1) t(end)]); grid on;

figure()
plot(t, filter_coil);
title('Coil Output for LVDT After Demodulation & Low Pass Filter');
xlabel('Time [s]'); ylabel('Voltage [V]');
xlim([t(1) t(end)]); ylim([0 2]); grid on;


% 2.2.4 Deflect the beam a small amount and release it capturing the step
% response

figure()
plot(tResponse, coil_response);
title('LVDT Deflection of a Beam Step Response');
xlabel('Time [s]'); ylabel('Voltage [V]');
xlim([tResponse(1) tResponse(end)]); grid on;

% 2.2 a.)

weight = xlsread('Calibration of LVDT.xlsx', 'A2:A12');
increasing_V = xlsread('Calibration of LVDT.xlsx', 'B2:B12');
decreasing_V = xlsread('Calibration of LVDT.xlsx', 'C2:C12');

distance = xlsread('Calibration of LVDT.xlsx', 'E2:E14');
lvdt_volt = xlsread('Calibration of LVDT.xlsx', 'F2:F14');

figure()
plot(weight, increasing_V, '*-', weight, decreasing_V, '*-');
title('Calibration of LVDT using Weights')
xlabel('Weight [g]'); ylabel('Voltage [mV]');
grid on; legend('Increasing', 'Decreasing')

figure()
plot(distance, lvdt_volt, '*-'); 
title('Calibration of LVDT using Micrometer')
xlabel('Distance [mili inch]'); ylabel('Voltage [mV]');
grid on; 

w = polyfit(weight, increasing_V, 1);
dist = polyfit(distance, lvdt_volt, 1);

% b.) 
weight_Sen     = w(1);              % mV/g
dist_Sen       = dist(1)*1000;           % mV/in
inchPerGram    = weight_Sen/dist_Sen;
dist_trav      = 13.3*inchPerGram;  % weight of the bucket
lbm            = 13.3/453.592;      % lbm 
lbf            = lbm*32.174;        %
beamSpring     = lbf/dist_trav;   

% c.) 
dist_Sen = dist_Sen/1000;
disp(strcat({'The system sensitivity in V/in is '}, num2str(dist_Sen, 4), ' V/in'))

% d.) 

th = 0.05;
[peakLoc, peakMag] = peakfinder(coil_response, th);
peakLoc(1)         = [];
peakMag(1)         = [];

n_peaks                 = length(peakLoc);
lvdt_Td                 = 1/((tResponse(peakLoc(50)) - tResponse(peakLoc(1)))/49);
lvdt_wd                 = lvdt_Td*2*pi;                        % wd

figure()
plot(tResponse, coil_response, tResponse(peakLoc), peakMag, 'd')
xlabel('Time [s]')
ylabel('LVDT')
title('Frequency Peaks of LVDT Voltage vs. Time')

for ii = 1:length(peakLoc)
    y(ii)          = log(peakMag(1)/peakMag(ii));
end

n = 0:length(y)-1;
damp_ratio         = zeros(1, length(peakLoc));

for jj = 1:length(peakLoc)
    num            = ((1/length(n))*log(peakMag(1)/peakMag(jj)));
    damp_ratio(jj) = num/(sqrt(4*pi^2 + num^2));
end

lvdt_zeta               = mean(damp_ratio);
lvdt_wn                 = lvdt_wd/(sqrt(1 - lvdt_zeta^2));    % wd = wn*sqrt(1-zeta^2)

% e.) effective mass of the system
meff     = beamSpring/lvdt_wn^2; %lbfs2/in

% part f
R = 100000
C=.047*10^-6
sys = tf([1],[ R*C 1])
bode(sys)

%% LVDT Frequency Response 
close all; clear all;

freq_amp   = xlsread('LVDT Frequency Response', 'Amplitude', 'A10:A50010');
amplitude  = xlsread('LVDT Frequency Response', 'Amplitude', 'B10:B50010');


freq_phase = xlsread('LVDT Frequency Response', 'Phase', 'A10:A50010');
phase      = xlsread('LVDT Frequency Response', 'Phase', 'B10:B50010');

figure()
subplot(2,1,1)
semilogx(freq_amp*2*pi, amplitude);
ylabel('Amplitude [db]')
xlim([freq_amp(1) freq_amp(end)*2*pi])
title('Theoretical Bode Plot')

subplot(2,1,2)
semilogx(freq_phase*2*pi, phase)
ylabel('Phase [\circ]')
xlabel('Frequency [rad/s]')
ylim([-400 400])
xlim([freq_phase(1) freq_phase(end)*2*pi])


wb1 = 295*2*pi;              % hz
wb2 = 22.07e3*2*pi;          % hz
Rp = 408;               % ohms
Rs = 162;               % ohms
Rm = 1e6;               % ohms

Lp = Rp/wb1;
L0 = (Rm + 2*Rs)/(wb2*2);

% c.) 

Ks = 10^(-22.57/20);          % db
x = 0.1;

Km = (Ks*Rp*(Rm + 2*Rs))/(2*Rm*x);

num = [((2*Rm*Km*x)/(Rp*(Rm + 2*Rs))) 0];
den = [( (Lp/Rp) * (2*L0/(Rm + 2*Rs)) ) (Lp/Rp + (2*L0/(Rm + 2*Rs))) 1];
sys = tf(num, den);

figure()
bode(sys)






