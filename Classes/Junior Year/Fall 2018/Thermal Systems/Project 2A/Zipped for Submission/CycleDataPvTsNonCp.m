
function [Pressure_Cycle,SpecificVolume_Cycle,Temperature_Cycle,SpecificEntropy_Cycle] = CycleDataPvTsNonCp(PR,MaxT)

% Thermal Systems - Project 2A   
% Function to obtain P-v and T-s data for non constant cp

%%%%%%
% Definition of Constants
%%%%%%
Air_Data = xlsread('air_data1.xls');
Temperatures    = Air_Data(:,1);
Specific_Enthalpy = Air_Data(:,3);
Specific_Entropy = Air_Data(:,4);

R=0.287;     % Constant

%%%%%%
% State Calculations
%%%%%%

% State 1 % Inlet of Compressor 
P1 = 100;    % Constant
T1 = 300;    % Constant
v1 = (R*T1)/P1; 
s1 = interp1(Temperatures,Specific_Entropy,T1);
st1 = interp1(Temperatures,Specific_Entropy,T1);
h1 = interp1(Specific_Entropy,Specific_Enthalpy,s1);

% State 2 % Inlet of Combuster  
P2 = PR*P1; 
s2 = s1;
st2 = st1 + R*log(PR);
T2 = interp1(Specific_Entropy,Temperatures,st2);
h2 = interp1(Temperatures,Specific_Enthalpy,T2);
v2 = R*T2/P2;

% State 3 % Inlet of Turbine
P3 = P2; % Isobaric Process
s3 = interp1(Temperatures,Specific_Entropy,MaxT)-st2+s2;
st3 = interp1(Temperatures,Specific_Entropy,MaxT);
v3 = (R * MaxT)/P3; 
h3 = interp1(Temperatures,Specific_Enthalpy,MaxT);


% State 4 % Inlet of Nozzle
s4 = s3;
h4 = h3 - (h2-h1);
T4 = interp1(Specific_Enthalpy,Temperatures,h4);
st4 = interp1(Temperatures,Specific_Entropy,T4);
P4 = P3*(exp((st4-st3)/R));
v4 = R*T4/P4;

% State 5 % Outlet of Nozzle
s5 = s4;
P5 = P1; 
st5 = st4 + R*log(P5/P4);
T5 = interp1(Specific_Entropy,Temperatures,st5) ;
v5 = (R*T5)/P5; 
h5 = interp1(Temperatures,Specific_Enthalpy,T5);


%%%%%%
% Process Calculations
%%%%%%

% Compressor - Isentropic
% Combustor  - Isobaric
% Turbine    - Isentropic
% Nozzle     - Isentropic


% Compressor Process - Isentropic - P as a function of v

P_Compressor = linspace(P1,P2,1000)';
s_Compressor = zeros(1000,1);
T_Compressor = linspace(T1,T2,1000)';
v_Compressor = zeros(1000,1);

for index = 1:1000
    s_Compressor(index) = s1;
    v_Compressor(index) = R*T_Compressor(index)/P_Compressor(index);
end

% Combustor Process - Isobaric - T as a function of s

P_Combustor = zeros(1000,1);
v_Combustor = linspace(v2,v3,1000)';
T_Combustor = linspace(T2,MaxT,1000)';
s_Combustor = zeros(1000,1);

for index = 1:1000
    P_Combustor(index) = P2;
    s_Combustor(index) = (interp1(Temperatures,Specific_Entropy,T_Combustor(index)))-st2+s2;
end

% Turbine Process - Isentropic - P as a function of v
s_Turbine = linspace(s3,s4,1000)';
v_Turbine = zeros(1000,1);
T_Turbine = linspace(MaxT,T4,1000)';
P_Turbine = linspace(P3,P4,1000)';


for index = 1:1000
    v_Turbine(index) = R*T_Turbine(index)/P_Turbine(index);
end


% Nozzle Process - Isentropic - P as a function of v

s_Nozzle = zeros(1000,1);
v_Nozzle = zeros(1000,1);
T_Nozzle = linspace(T4,T5,1000)';
P_Nozzle = linspace(P4,P5,1000)';

for index = 1:1000
    s_Nozzle(index) = s4;
    v_Nozzle(index) = R*T4/P_Nozzle(index);
end

% Reset Process - 4 to 1 - Isobaric

v_Reset = linspace(v5,v1,1000)';
P_Reset = zeros(1000,1);
T_Reset = linspace(T5,T1,1000)';
s_Reset = zeros(1000,1);

for index = 1:1000
    P_Reset(index) = P5;
    s_Reset(index) = interp1(Temperatures,Specific_Entropy,T_Reset(index))-st1+s1;
end

% Data Compiled for Output - 1000 points per process
Pressure_Cycle        = [P_Compressor;P_Combustor;P_Turbine;P_Nozzle;P_Reset];
SpecificVolume_Cycle  = [v_Compressor;v_Combustor;v_Turbine;v_Nozzle;v_Reset];
Temperature_Cycle     = [T_Compressor;T_Combustor;T_Turbine;T_Nozzle;T_Reset];
SpecificEntropy_Cycle = [s_Compressor;s_Combustor;s_Turbine;s_Nozzle;s_Reset];
end















