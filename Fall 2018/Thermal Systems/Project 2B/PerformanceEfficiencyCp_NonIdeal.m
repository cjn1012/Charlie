
clear all
close all

% Thermal Systems - Project 2A   
% Function to obtain performance and efficient of a jet engine with particular parameters for constant cp


T1 = 300
P1 = 100
PR = 30
MaxT = 1500
%%%%%%
% Definition of Constants
%%%%%%

R=0.287;     % Constant
n=1.4;       % Constant
k = (1/0.4); % Constant
Cp = 1.005; % kJ/Kg

%%%%%%
% State Calculations for Ideal
%%%%%%

% State 1 % Inlet of Compressor 
v1 = (R*T1)/P1; 
s1 = 1;
C=(P1*(v1^n));  % Isentropic Constant
 
% State 2 % Inlet of Combuster  
P2 = PR*P1; 
v2 = (C/P2)^(1/n);   
T2= (P2*v2)/R;  
s2 = 1;

% State 3 % Inlet of Turbine   
P3 = P2; % Isobaric Process
v3 = (R * MaxT)/P3; 
s3 = s2 + Cp*log(MaxT/T2);

% State 4 % Inlet of Nozzle
T4 = MaxT - (T2-T1);   
v4 = ((P3 * (v3^n))/(R*T4))^(k); 
P4 = (R * T4)/v4;
s4 = s3;

% State 5 % Outlet of Nozzle
P5 = P1; 
v5 = ((P4*(v4^n))/P5)^(1/n); 
T5 = (P5*v5)/R;   
s5 = s4;


%%%%%%
% State Calculations for Non-Ideal
%%%%%%
T1n = T1;
P1n = P1;
MaxTn = MaxT;
PRn = PR;

Eff_Comp = .95;
Eff_Turb = .9;
Eff_Nozz = .9;

% State 1 % Inlet of Compressor 
v1n = (R*T1n)/P1n; 
s1n = 1;

 
% State 2 % Inlet of Combuster  
P2n = PRn*P1n;   
T2n= (T2-T1n)/Eff_Comp + T1n;  
s2n = s1n + Cp*log(T2n/T1n)-R*log(P2n/P1n);

% State 3 % Inlet of Turbine   
P3n = P2n; % Isobaric Process
v3n = (R * MaxTn)/P3n; 
s3n = s2n + Cp*log(MaxTn/T2n);

% State 4 % Inlet of Nozzle
T4n = -Eff_Turb * (MaxTn-T4) + MaxTn ; 
s4n = s3n + Cp*log(T4n/MaxTn)

% State 5 % Outlet of Nozzle
P5n = P1n; 
T5n = -Eff_Nozz * (T4n-T5) + T4n ;   
s5n = s4n + Cp*log(T5n/T4n);






%%%%%
% Calculating the variables needed to find the 
% overall efficiency and performance of the jet engine
%%%%%


Q_Add    = Cp*(MaxT-T2);
Q_Remove = Cp*(T5-T1);

W_Compressor = Cp*(T2-T1);
W_Turbine    = Cp*(T4-MaxT);


Performance = W_Compressor + Q_Add - Q_Remove + W_Turbine;
Thermal_Efficiency = Performance/Q_Add; % Dimensionless Parameter


linesy = [T1,T2,MaxT,T4,T5,T1];
linesx = [s1,s2,s3,s4,s5,s1];
linesnx = [s1n,s2n,s3n,s4n];
linesny = [T1n,T2n,MaxTn,T4n];

plot(linesx,linesy)
hold on
plot(linesnx,linesny)
















