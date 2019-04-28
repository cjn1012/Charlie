clear all
close all

% Reading in Excel Load Cell File
MomentVoltage = xlsread('LoadCell_NoFairings.xlsx','Data' ,'A2:A50008'); % Second
Moment = MomentVoltage.*33405-.7733;

D1 = 16; % Length of pin to load cell
D2 = 43; % Length of pin to concentrated drag force

Force = (Moment*D1)/D2
Time = linspace(0,length(Force)/100,length(Force))

figure(1)
plot(Time,Force)




%% With Fairings
% Reading in Excel Load Cell File
MomentVoltage2 = xlsread('LoadCell_Fairings.xlsx','Data' ,'A2:A50008'); % Second
Moment2 = MomentVoltage2.*33405-.7733;

D1 = 16; % Length of pin to load cell
D2 = 43; % Length of pin to concentrated drag force

Force2 = (Moment2*D1)/D2;
Time2 = linspace(0,length(Force2)/100,length(Force2));

figure(2)
plot(Time2,Force2)






