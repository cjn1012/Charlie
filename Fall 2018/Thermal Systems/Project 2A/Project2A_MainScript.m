% Script to create all the graphs needed for report using
% functions within this folder
clear all;
close all;


% Calculating Performance and Efficiency for Varying PR

T1   = 300;
P1   = 100;
MaxT1 = 1500;
PR   = linspace(1,50,50)';

EfficiencyPRNonCp  = zeros(50,1);
PerformancePRNonCp = zeros(50,1);
EfficiencyPRCp  = zeros(50,1);
PerformancePRCp = zeros(50,1);
for i = 1:50
    [PerformancePRNonCp(i),EfficiencyPRNonCp(i)] = PerformanceEfficiencyNonCp(T1,P1,PR(i),MaxT1);
    [PerformancePRCp(i),EfficiencyPRCp(i)] = PerformanceEfficiencyCp(T1,P1,PR(i),MaxT1);
end



% Calculating Performance and Efficiency for Varying MaxT

T1   = 300;
P1   = 100;
PR2   = 30;
MaxT = linspace(1000,2000,50)';

EfficiencyMaxTNonCp  = zeros(50,1);
PerformanceMaxTNonCp = zeros(50,1);
EfficiencyMaxTCp  = zeros(50,1);
PerformanceMaxTCp = zeros(50,1);
for i = 1:50
    [PerformanceMaxTNonCp(i),EfficiencyMaxTNonCp(i)] = PerformanceEfficiencyNonCp(T1,P1,PR2,MaxT(i));
    [PerformanceMaxTCp(i),EfficiencyMaxTCp(i)] = PerformanceEfficiencyCp(T1,P1,PR2,MaxT(i));
end


% Calculating P,v,T and s data points for various parameters

% Data for the baseline - average PR and MaxT
PR3   = 25;
MaxT2 = 1500;
[PressuresBase,SpecificVolumesBase,TemperaturesBase,SpecificEntropysBase] = CycleDataPvTs(PR3,MaxT2);

% Data for increased MaxT - average PR
PR4   =25;
MaxT3 = 2000;
[PressuresMaxT,SpecificVolumesMaxT,TemperaturesMaxT,SpecificEntropysMaxT] = CycleDataPvTs(PR4,MaxT3);

% Data for increased PR - average MaxT
PR5   = 35;
MaxT4 = 1500;
[PressuresPR,SpecificVolumesPR,TemperaturesPR,SpecificEntropysPR] = CycleDataPvTs(PR5,MaxT4);



% Bonus Problem, calculating efficiency and performance as a function of PR and MaxT

[PRx,MaxTy] = meshgrid(PR,MaxT);

Performance = zeros(50,50);
Thermal_Efficiency = zeros(50,50);

for i = 1:50
    for j = 1:50
        [Performance(i,j),Thermal_Efficiency(i,j)] = PerformanceEfficiencyNonCp(T1,P1,PR(i),MaxT(j));
    end
end


%%%%%
% Plotting Data
%%%%%

figure(1)

plot(PR,EfficiencyPRNonCp,'r')
hold on
plot(PR,EfficiencyPRCp,'b')
ylabel('Efficiency ([])','FontSize',20)
set(gca,'fontsize',18)
xlabel('Compressor Pressure Ratio ([])','FontSize',20)
set(gca,'fontsize',18)
lgd = legend('\color{red} Non-Constant Cp','\color{blue} Constant Cp','Location','northwest');
lgd.FontSize = 14;
hold off



figure(2)
plot(PR,PerformancePRNonCp,'r')
hold on
plot(PR,PerformancePRCp,'b')
ylabel('Performance ([])','FontSize',20)
set(gca,'fontsize',18)
xlabel('Compressor Pressure Ratio ([])','FontSize',20)
set(gca,'fontsize',18)
lgd = legend('\color{red} Non-Constant Cp','\color{blue} Constant Cp','Location','northwest');
lgd.FontSize = 14;
hold off

figure(3)
plot(MaxT,EfficiencyMaxTNonCp,'r')
hold on
plot(MaxT,EfficiencyMaxTCp,'b')
ylabel('Efficiency ([])','FontSize',20)
set(gca,'fontsize',18)
xlabel('Combustion Maximum Temperature [K])','FontSize',20)
set(gca,'fontsize',18)
lgd = legend('\color{red} Non-Constant Cp','\color{blue} Constant Cp','Location','northwest');
lgd.FontSize = 14;
hold off

figure(4)
plot(MaxT,PerformanceMaxTNonCp,'r')
hold on
plot(MaxT,PerformanceMaxTCp,'b')
ylabel('Efficiency ([])','FontSize',20)
set(gca,'fontsize',18)
xlabel('Combustion Maximum Temperature [K])','FontSize',20)
set(gca,'fontsize',18)
lgd = legend('\color{red} Non-Constant Cp','\color{blue} Constant Cp','Location','northwest');
lgd.FontSize = 14;
hold off

figure(5)
mesh(PRx,MaxTy,Thermal_Efficiency)

figure(6)
mesh(PRx,MaxTy,Performance)


figure(7)

plot(SpecificVolumesBase,PressuresBase,'k',SpecificVolumesMaxT,PressuresMaxT,'r',SpecificVolumesPR,PressuresPR,'b')
ylabel('Pressure [KPa])','FontSize',20)
set(gca,'fontsize',18)
xlabel('Specific Volume [m^3/kg])', 'Interpreter', 'none','FontSize',20)
set(gca,'fontsize',18)
lgd = legend('\color{black} Baseline Cycle','\color{red} Increased Max Temperature Cycle','\color{blue} Increased Compressor Pressure Ratio Cycle');
lgd.FontSize = 12;
hold off

figure(8)

plot(SpecificEntropysBase,TemperaturesBase,'k',SpecificEntropysMaxT,TemperaturesMaxT,'r',SpecificEntropysPR,TemperaturesPR,'b')
ylabel('Temperature [K])','FontSize',20)
set(gca,'fontsize',18)
xlabel('Specific Entropy [KJ/KG])','FontSize',20)
set(gca,'fontsize',18)
lgd = legend('\color{black} Baseline Cycle','\color{red} Increased Max Temperature Cycle','\color{blue} Increased Compressor Pressure Ratio Cycle');
lgd.FontSize = 12;
hold off












