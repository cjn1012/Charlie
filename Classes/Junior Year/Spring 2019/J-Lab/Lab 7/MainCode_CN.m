clear all
close all

% Reading in Excel Load Cell File
MomentVoltage = xlsread('LoadCell_NoFairings.xlsx','Data' ,'A2:A50008'); % Second
Moment = MomentVoltage.*33405-.7733;

D1 = 16; % Length of pin to load cell
D2 = 43; % Length of pin to concentrated drag force

Force = (Moment*D1)/D2 ;
Force_Cal = 0 - Force(1);
Force = Force + Force_Cal;
Force = Force/.225;
Time = linspace(0,length(Force)/100,length(Force))-2;

figure(1)
plot(Time,Force)
xlabel('Time (s)')
ylabel('Drag Force (N)')
title('Drag Force - No Fairings')
xlim([0,14])



%% With Fairings
% Reading in Excel Load Cell File
MomentVoltage2 = xlsread('LoadCell_Fairings.xlsx','Data' ,'A2:A50008'); % Second
Moment2 = MomentVoltage2.*33405-.7733;

D1 = 16; % Length of pin to load cell
D2 = 43; % Length of pin to concentrated drag force

Force2 = (Moment2*D1)/D2;
Force2_Cal = 0 - Force2(3);
Force2 = Force2 + Force2_Cal;
Force2 = Force2/.225
Time2 = linspace(0,length(Force2)/100,length(Force2))-5;

figure(2)
plot(Time2,Force2)
xlabel('Time (s)')
ylabel('Drag Force (N)')
title('Drag Force - With Fairings')
xlim([0,14])




figure(3)
subplot(2,1,1)
plot(Time,Force)
xlabel('Time (s)')
ylabel('Drag Force (N)')
title('Drag Force - No Fairings')
xlim([0,14])

subplot(2,1,2)
plot(Time2,Force2)
xlabel('Time (s)')
ylabel('Drag Force (N)')
title('Drag Force - With Fairings')
xlim([0,14])





