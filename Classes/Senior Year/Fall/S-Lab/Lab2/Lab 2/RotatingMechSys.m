clear all; close all; clc; 

%% Position and Velocity Measurements

potSens = .015;                                      % V/degree   
tacSens = 7/1000;                                    % mv/rpm

angPos  = xlsread('PosVelMes.xlsx', 'B8:B100007');   % Pot
angPos  = (angPos/potSens)*((2*pi)/360);             % rad
angVel  = xlsread('PosVelMes.xlsx', 'A8:A100007');   % tach
angVel  = (angVel/tacSens)*2*pi*(1/60);              % rad/s

t       = 0:1.201e-5:(1.201e-5*1e5);
t       = t(1:end-1);

% a.)i. 
f1=figure(1)
yyaxis left 
plot(t, angPos);
xlabel('Time [s]'); ylabel('Angular Position [rad]');

yyaxis right 
plot(t, angVel);
ylabel('Angular Velocity [rad/s]');
legend('Potentiometer', 'Tachometer');
title('Angular Position and Velocity vs Time')

%a.)ii.
f2=figure(2)
plot(angVel, angPos);
xlabel('Potentiometer [rad]'); ylabel('Tachometer [rad/s]');
title('Angular Velocity vs. Angular Position');

%d.) 
intPos = cumtrapz(t, angVel);
f3=figure(3)
plot(t, angPos - mean(angPos(1:5000)), t, intPos);
xlabel('Time [s]'); ylabel('Angular Position [rad]');
title('Integrated & Measured Angular Position');
legend('Measured \theta', 'Integrated \theta');

%% Second-Order Parameter Identification

Jtach              = 1.32e-4*(1/16);
Jpot               = 2.44e-3*(1/16);
convert            = 1/(32.174*12);
DenSteel           = 0.285;                         % lbm/in^3

% b.)
A1                 = ((pi*.065^2)/4)/(12^2);        % ft^2
L1                 = 5/12;                          % ft
L2                 = 4.75/12;                       % ft
K1                 = 0.253;                         % lbf-ft/rad
E1                 = (K1*L1)/A1;                    % lbf/rad
Kwire              = ((E1*A1)/L2);                  % lbf-ft/rad
Kwire              = ((E1*A1)/L2)*12;               % lbf-in/rad
 
WireVol            = pi*L2*(.065/2)^2;              % in^3
WireMass           = DenSteel*WireVol;              % lbm
RodVol             = pi*(5.6)*(.25/2)^2;            % in^3
RodMass            = DenSteel*RodVol;               % lbm

% c.)
% Smoothing Angular Position data and finding peaks in data
th                 = .05;
smooth             = 9;
[angPos, smooth]   = wsmooth(angPos, t, smooth);
[peakLoc, peakMag] = peakfinder(angPos - mean(angPos(1:5000)), th);
peakLoc(1)         = [];
peakMag(1)         = [];

% Damped Natural Frequency
n_peaks            = length(peakLoc);
Td                 = n_peaks/(t(peakLoc(end)) - t(peakLoc(1)));
AngPos_wd          = Td*2*pi;                        % wd

f4=figure(4)
plot(t, angPos - mean(angPos(1:5000)), t(peakLoc), peakMag, 'd')
xlabel('Time [s]')
ylabel('Angular Position [\theta]')
title('Frequency Peaks of Angular Position vs. Time')

for ii = 1:length(peakLoc)
    y(ii)          = log(peakMag(1)/peakMag(ii));
end

n = 0:length(y)-1;
damp_ratio         = zeros(1, length(peakLoc));

for jj = 1:length(peakLoc)
    num            = ((1/length(n))*log(peakMag(1)/peakMag(jj)));
    damp_ratio(jj) = num/(sqrt(4*pi^2 + num^2));
end

AngPos_zeta        = mean(damp_ratio);
AngPos_wn          = AngPos_wd/(sqrt(1 - AngPos_zeta^2));    % wd = wn*sqrt(1-zeta^2)


% Smoothing Angular Velocity data and finding peaks in data
th                 = .015;
smooth             = 7;
[angVel, smooth]   = wsmooth(angVel, t, smooth);
[peakLoc, peakMag] = peakfinder(angVel, th);
peakLoc(1)         = [];
peakMag(1)         = [];

% Damped Natural Frequency
n_peaks            = length(peakLoc);
Td                 = n_peaks/(t(peakLoc(end)) - t(peakLoc(1)));
AngVel_wd          = Td*2*pi;                        % wd

f5=figure(5)
plot(t, angVel, t(peakLoc), peakMag, 'd')
xlabel('Time [s]')
ylabel('Angular Velocity [\theta/s]')
title('Frequency Peaks of Angular Velocity vs. Time')

for ii = 1:length(peakLoc)
    y(ii)          = log(peakMag(1)/peakMag(ii));
end

n = 0:length(y)-1;
damp_ratio         = zeros(1, length(peakLoc));

for jj = 1:length(peakLoc)
    num            = ((1/length(n))*log(peakMag(1)/peakMag(jj)));
    damp_ratio(jj) = num/(sqrt(4*pi^2 + num^2));
end

AngVel_zeta        = mean(damp_ratio);
AngVel_wn          = AngVel_wd/(sqrt(1 - AngVel_zeta^2));    % wd = wn*sqrt(1-zeta^2)

% d.)
J_Exp              = ((Kwire)/((AngVel_wn)^2));
B                  = (2*AngVel_zeta*Kwire)/AngVel_wn;


% e.) 
Gsteel             = 16.562e8 / (12^2);                      % lbf/in^2
J_gear             = convert*(1/2)*(.268)*(2.5/2)^2;
J_C1               = convert*(1/2)*(.0249)*(.744/2)^2;
J_C2               = convert*(1/2)*(.0253)*(.75/2)^2;
J_Wire             = convert*(1/2)*(WireMass)*(.065/2)^2;
J_Rod              = convert*(1/2)*(RodMass)*(.25/2)^2;
J_theor            = J_gear + Jpot + Jtach + J_C1 + J_C2 + J_Wire + J_Rod;    % lbf-in-s^2/rad

Ktheor             = Gsteel*(J_Wire/(L2*12));

%% Integration of Velocity and Signal Drift

potPos   = xlsread('IntegratorWithDrift.xlsx', 'B8:B100007');   % Pot
tachVel  = xlsread('IntegratorWithDrift.xlsx', 'A8:A100007');   % tach

R                  = 10e3;
C                  = 10e-6;
tacSens2 = 0.007*60/(2*pi);
potPos  = (angPos/potSens)*((2*pi)/360);             % rad
tachVel = (tachVel/tacSens2)*(R*C); % Not really sure if i need this part: *(2*pi)*(1/60)

f6=figure(6)
plot(t, tachVel - tachVel(1), t, potPos - potPos(1))
xlabel('Time [s]'); ylabel('Angular Position [rad]');
legend('Integrator Op-Amp Signal', 'Potentiometer');
title('Angle of Potentiometer and Tachometer with Drift vs. Time')

%% Integration of Velocity Without Signal Drift

pot     = xlsread('IntegratorWithoutDrift.csv', 'B8:B100007');   
tac     = xlsread('IntegratorWithoutDrift.csv', 'A8:A100007');

Rf      = 1.01e6; 
Ri      = 9.97e3;
Cf      = 10e-6;

num     = -Rf/Ri;
den     = [Rf*Cf 1];
sys     = tf(num, den);

in =[.802 .803 .8];
out =[8.8 .962 .120];
dt =[1.52 .161 .016];
Freq =[.15 1.5 15]*2*pi;
period =[6.66 .667 .0668];

for i=1:length(dt)
    mag(i)=20*log10(out(i)/in(i));
    phase(i)=(dt(i)/period(i))*360;
end


f7=figure(7)
hold on

semilogx(Freq,mag,'k--')
axis([10e-4 10e4 -10 45])

bode(sys)

semilogx(Freq/(2*pi),phase,'k--')
axis([10e-3 10e2 80 180])

legend('Theoretical','Experimental Data')

grid on





% b.) 

pottheta  = (pot/.015)*(pi/180);             % rad
tactheta = (tac/.007)*(Ri*Cf)*2*pi*(1/60); % Not really sure if i need this part: *(2*pi)*(1/60)
pottheta  = pottheta - pottheta(1);   % zeroing
tactheta  = tactheta - tactheta(1);

f8=figure(8)
plot(t,pottheta,t,tactheta)
legend('Potentiometer Angle Disp','Tachometer Angle Disp')
xlabel('Time (seconds)')
ylabel('Angle Displacement (rad)')
title('Angle of Potentiometer and Tachometer without Drift vs. Time')
grid on


%C

O=9.7;
I=.125;
BP_Mag=20*log10(O/I);
fit = polyfit(log10(Freq/(2*pi)),mag,1);
m1=fit(1);
b1=fit(2);
BPfreq=10^((BP_Mag-b1)/m1)

tau_theo = Rf*Cf
gainvalue = Rf/Ri
Tau_exper = (1/BPfreq)*(1/(2*pi))
gainvalue_T = 10^(BP_Mag/20)


% Save all figures
saveas(f1,'f1','png')
saveas(f2,'f2','png')
saveas(f3,'f3','png')
saveas(f4,'f4','png')
saveas(f5,'f5','png')
saveas(f6,'f6','png')
saveas(f7,'f7','png')
saveas(f8,'f8','png')





