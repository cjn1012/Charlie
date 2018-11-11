% Script to create all the graphs needed for report using
% functions within this folder
clear all;
close all;


%%%%%%%%%%
% Calculating Performance and Efficiency for Varying PR
%%%%%%%%%%

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

%%%%%%%%%%
% Calculating Performance and Efficiency for Varying MaxT
%%%%%%%%%%

T1   = 300;
P1   = 100;
PR2   = 25;
MaxT = linspace(1000,2000,50)';

EfficiencyMaxTNonCp  = zeros(50,1);
PerformanceMaxTNonCp = zeros(50,1);
EfficiencyMaxTCp  = zeros(50,1);
PerformanceMaxTCp = zeros(50,1);
for i = 1:50
    [PerformanceMaxTNonCp(i),EfficiencyMaxTNonCp(i)] = PerformanceEfficiencyNonCp(T1,P1,PR2,MaxT(i));
    [PerformanceMaxTCp(i),EfficiencyMaxTCp(i)] = PerformanceEfficiencyCp(T1,P1,PR2,MaxT(i));
end

%%%%%%%%%%
% Calculating P,v,T and s data points for various parameters constant cp
%%%%%%%%%%

% Data for the baseline - average PR and MaxT
PR3   = 25;
MaxT2 = 1500;
[PressuresBase,SpecificVolumesBase,TemperaturesBase,SpecificEntropysBase] = CycleDataPvTsCp(PR3,MaxT2);

% Data for increased MaxT - average PR
PR4   = 25;
MaxT3 = 2000;
[PressuresMaxT,SpecificVolumesMaxT,TemperaturesMaxT,SpecificEntropysMaxT] = CycleDataPvTsCp(PR4,MaxT3);

% Data for increased PR - average MaxT
PR5   = 40;
MaxT4 = 1500;
[PressuresPR,SpecificVolumesPR,TemperaturesPR,SpecificEntropysPR] = CycleDataPvTsCp(PR5,MaxT4);


%%%%%%%%%%
% Calculating P,v,T and s data points for various parameters non constant cp
%%%%%%%%%%

% Data for the baseline - average PR and MaxT
PR3   = 25;
MaxT2 = 1500;
[PressuresBase2,SpecificVolumesBase2,TemperaturesBase2,SpecificEntropysBase2] = CycleDataPvTsNonCp(PR3,MaxT2);

% Data for increased MaxT - average PR
PR4   = 25;
MaxT3 = 2000;
[PressuresMaxT2,SpecificVolumesMaxT2,TemperaturesMaxT2,SpecificEntropysMaxT2] = CycleDataPvTsNonCp(PR4,MaxT3);

% Data for increased PR - average MaxT
PR5   = 40;
MaxT4 = 1500;
[PressuresPR2,SpecificVolumesPR2,TemperaturesPR2,SpecificEntropysPR2] = CycleDataPvTsNonCp(PR5,MaxT4);


% Bonus Problem, calculating efficiency and performance as a function of PR and MaxT

Res = 10;

PR2   = linspace(1,50,Res)';
MaxT3 = linspace(1000,2000,Res)';
T1 = 300;
P1 = 100;
[PRx,MaxTy] = meshgrid(PR2,MaxT3);
Performance = zeros(Res,Res);
Thermal_Efficiency = zeros(Res,Res);

for i = 1:Res
    for j = 1:Res
        [Performance(i,j),Thermal_Efficiency(i,j)] = PerformanceEfficiencyNonCp(T1,P1,PR2(i),MaxT3(j));
    end
end

%%%%%%%%%%
% Plotting Data for Efficiency and Performance
%%%%%%%%%%

% Plot of Efficiency vs. Pressure Ratio
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


% Plot of Performance vs. Pressure Ratio
figure(2)
plot(PR,PerformancePRNonCp,'r')
hold on
plot(PR,PerformancePRCp,'b')
ylabel('Performance ([m^2/s^2])','FontSize',20)
set(gca,'fontsize',18)
xlabel('Compressor Pressure Ratio ([])','FontSize',20)
set(gca,'fontsize',18)
lgd = legend('\color{red} Non-Constant Cp','\color{blue} Constant Cp','Location','northwest');
lgd.FontSize = 14;
hold off

% Plot of Efficiency vs. Maximum Temperature
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

% Plot of Performance vs. Maximum Temperature
figure(4)
plot(MaxT,PerformanceMaxTNonCp,'r')
hold on
plot(MaxT,PerformanceMaxTCp,'b')
ylabel('Performance ([m^2/s^2])','FontSize',20)
set(gca,'fontsize',18)
xlabel('Combustion Maximum Temperature [K])','FontSize',20)
set(gca,'fontsize',18)
lgd = legend('\color{red} Non-Constant Cp','\color{blue} Constant Cp','Location','northwest');
lgd.FontSize = 14;
hold off


%%%%%%%%%%
% Bonus 3D Plots
%%%%%%%%%%

figure(5)
mesh(PRx,MaxTy,Thermal_Efficiency)
xlabel('Pressure Ratio [ ])','Rotation',15,'FontSize',15)
set(gca,'fontsize',14)
ylabel('Maximum Temperature [K])', 'Rotation',-25,'FontSize',15)
set(gca,'fontsize',14)
zlabel('Efficiency [ ])', 'Interpreter', 'none','FontSize',15)
set(gca,'fontsize',14)
hold off

figure(6)
mesh(PRx,MaxTy,Performance)
xlabel('Pressure Ratio [ ])','Rotation',15,'FontSize',15)
set(gca,'fontsize',14)
ylabel('Maximum Temperature [K])', 'Rotation',-25,'FontSize',15)
set(gca,'fontsize',14)
zlabel('Performance [m^2/s^2])','FontSize',15)
set(gca,'fontsize',14)
hold off


%%%%%%%%%%
% Cycle Graphs Constant Cp
%%%%%%%%%%

% Plot of Pressure vs. Specific Volume
figure(7)
plot(SpecificVolumesBase,PressuresBase,'k',SpecificVolumesMaxT,PressuresMaxT,'r',SpecificVolumesPR,PressuresPR,'b')
ylabel('Pressure [KPa])','FontSize',20)
set(gca,'fontsize',18)
xlabel('Specific Volume [m^3/kg])','FontSize',20)
set(gca,'fontsize',18)
lgd = legend('\color{black} Baseline Cycle','\color{red} Increased Max Temperature Cycle','\color{blue} Increased Compressor Pressure Ratio Cycle');
lgd.FontSize = 12;
hold off

% Plot of Temperature vs. Specific Entropy
figure(8)
plot(SpecificEntropysBase,TemperaturesBase,'k',SpecificEntropysMaxT,TemperaturesMaxT,'r',SpecificEntropysPR,TemperaturesPR,'b')
ylabel('Temperature [K])','FontSize',20)
set(gca,'fontsize',18)
xlabel('Specific Entropy [KJ/KG])','FontSize',20)
set(gca,'fontsize',18)
lgd = legend('\color{black} Baseline Cycle','\color{red} Increased Max Temperature Cycle','\color{blue} Increased Compressor Pressure Ratio Cycle');
lgd.FontSize = 12;
hold off


%%%%%%%%%%
% Cycle Graphs Non-Constant Cp
%%%%%%%%%%

% Plot of Pressure vs. Specific Volume
figure(9)

plot(SpecificVolumesBase2,PressuresBase2,'k',SpecificVolumesMaxT2,PressuresMaxT2,'r',SpecificVolumesPR2,PressuresPR2,'b')
ylabel('Pressure [KPa])','FontSize',20)
set(gca,'fontsize',18)
xlabel('Specific Volume [m^3/kg])','FontSize',20)
set(gca,'fontsize',18)
lgd = legend('\color{black} Baseline Cycle','\color{red} Increased Max Temperature Cycle','\color{blue} Increased Compressor Pressure Ratio Cycle');
lgd.FontSize = 12;
hold off

% Plot of Temperature vs. Specific Entropy
figure(10)
plot(SpecificEntropysBase2,TemperaturesBase2,'k',SpecificEntropysMaxT2,TemperaturesMaxT2,'r',SpecificEntropysPR2,TemperaturesPR2,'b')
ylabel('Temperature [K])','FontSize',20)
set(gca,'fontsize',18)
xlabel('Specific Entropy [KJ/kg])','FontSize',20)
set(gca,'fontsize',18)
lgd = legend('\color{black} Baseline Cycle','\color{red} Increased Max Temperature Cycle','\color{blue} Increased Compressor Pressure Ratio Cycle');
lgd.FontSize = 12;
hold off


% Plot of Temperature vs. Specific Entropy
figure(11)
plot(SpecificEntropysBase,TemperaturesBase,'k')
ylabel('Temperature [K])','FontSize',20)
set(gca,'fontsize',18)
xlabel('Specific Entropy [KJ/KG])','FontSize',20)
set(gca,'fontsize',18)
lgd = legend('\color{black} Baseline Cycle');
lgd.FontSize = 12;
hold off


