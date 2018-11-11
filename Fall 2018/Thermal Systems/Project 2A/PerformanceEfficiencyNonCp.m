
function [Performance,Thermal_Efficiency] = PerformanceEfficiencyNonCp(T1,P1,PR,MaxT)

% Thermal Systems - Project 2A   
% Function to obtain performance and efficient of a jet engine with particular parameters for non constant cp

%%%%%%
% Definition of Constants
%%%%%%
Air_Data = xlsread('air_data1.xls');
Temperatures    = Air_Data(:,1);
Internal_Energy = Air_Data(:,3);
Specific_Entropy = Air_Data(:,4);

R=0.287;     % Constant

%%%%%%
% State Calculations
%%%%%%

% State 1 % Inlet of Compressor 
v1 = (R*T1)/P1; 
s1 = interp1(Temperatures,Specific_Entropy,T1);
h1 = interp1(Specific_Entropy,Internal_Energy,s1);
st1 = s1;

% State 2 % Inlet of Combuster  
P2 = PR*P1; 
s2 = s1;
st2 = st1 + R*log(PR);
T2 = interp1(Specific_Entropy,Temperatures,st2);
h2 = interp1(Temperatures,Internal_Energy,T2);
v2 = R*T2/P2;

% State 3 % Inlet of Turbine
P3 = P2; % Isobaric Process
s3 = interp1(Temperatures,Specific_Entropy,MaxT)-st2+s2;
st3 = interp1(Temperatures,Specific_Entropy,MaxT);
v3 = (R * MaxT)/P3; 
h3 = interp1(Temperatures,Internal_Energy,MaxT);

% State 4 % Inlet of Nozzle
s4 = s3;
T4 = interp1(Specific_Entropy,Temperatures,s4);
st4 = interp1(Temperatures,Specific_Entropy,T4);
h4 = h3 - (h2-h1);
P4 = P3*(exp((st4-st3)/R));
v4 = R*T4/P4;

% State 5 % Outlet of Nozzle
s5 = s4;
P5 = P1; 
st5 = st4 + R*log(P5/P4);
T5 = interp1(Specific_Entropy,Temperatures,st5) ;
v5 = (R*T5)/P5; 
h5 = interp1(Temperatures,Internal_Energy,T5);


%%%%%
% Calculating the variables needed to find the 
% overall efficiency and performance of the jet engine
%%%%%

Q_Add    = h3 - h2;
Q_Remove = h5 - h1;

W_Compressor = h2 - h1;
W_Turbine    = h4 - h3;

Performance = W_Compressor + Q_Add - Q_Remove + W_Turbine;
Thermal_Efficiency = Performance/Q_Add;

end















