
function [Performance,Thermal_Efficiency] = PerformanceEfficiencyNonCp(T1,P1,PR,MaxT)

% Thermal Systems - Project 2A   
% Function to obtain performance and efficient of a jet engine with particular parameters

%%%%%%
% Definition of Constants
%%%%%%

R=0.287;     % Constant
n=1.4;       % Constant
k = (1/0.4); % Constant

%%%%%%
% State Calculations
%%%%%%

% State 1 % Inlet of Compressor 
v1 = (R*T1)/P1; 
C=(P1*(v1^n));  % Isentropic Constant
 
% State 2 % Inlet of Combuster  
P2 = PR*P1; 
v2 = (C/P2)^(1/n);   
T2= (P2*v2)/R;  

% State 3 % Inlet of Turbine   
P3 = P2; % Isobaric Process
v3 = (R * MaxT)/P3; 

% State 4 % Inlet of Nozzle
T4 = MaxT - (T2-T1);   
v4 = ((P3 * (v3^n))/(R*T4))^(k); 
P4 = (R * T4)/v4;

% State 5 % Outlet of Nozzle
P5 = P1; 
v5 = ((P4*(v4^n))/P5)^(1/n); 
T5 = (P5*v5)/R;   

%%%%%
% Calculating the variables needed to find the 
% overall efficiency and performance of the jet engine
%%%%%

Air_Data = xlsread('air_data1.xls');
Temperatures    = Air_Data(:,1);
Internal_Energy = Air_Data(:,2);

u_T1   = interp1(Temperatures,Internal_Energy,T1);
u_T2   = interp1(Temperatures,Internal_Energy,T2);
u_MaxT = interp1(Temperatures,Internal_Energy,MaxT);
u_T4   = interp1(Temperatures,Internal_Energy,T4);
u_T5   = interp1(Temperatures,Internal_Energy,T5);

Q_Add    = u_MaxT - u_T2;
Q_Remove = u_T5 - u_T1;

W_Compressor = u_T2 - u_T1;
W_Turbine    = u_T4 - u_MaxT;


Performance = W_Compressor + Q_Add - Q_Remove + W_Turbine;
Thermal_Efficiency = Performance/Q_Add; % Dimensionless Parameter

end















