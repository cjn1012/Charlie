
function [Performance,Thermal_Efficiency] = PerformanceEfficiencyCp(T1,P1,PR,MaxT)

% Thermal Systems - Project 2A   
% Function to obtain performance and efficient of a jet engine with particular parameters for constant cp

%%%%%%
% Definition of Constants
%%%%%%

R=0.287;     % Constant
n=1.4;       % Constant
k = (1/0.4); % Constant
Cp = 1.005; % kJ/Kg

%%%%%%
% State Calculations
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
s5 = s1;


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

end















