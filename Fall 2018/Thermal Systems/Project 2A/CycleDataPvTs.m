
function [Pressure_Cycle,SpecificVolume_Cycle,Temperature_Cycle,SpecificEntropy_Cycle] = CycleDataPvTs(PR,MaxT)

% Thermal Systems - Project 2A   
% Function to obtain P-v and T-s data

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
T1 = 300;       % Inlet Temperature
P1=100;         % Inlet Pressure
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
C2=(P3*(v3^n));

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


%%%%%%
% Process Calculations
%%%%%%

% Compressor - Isentropic
% Combustor  - Isobaric
% Turbine    - Isentropic
% Nozzle     - Isentropic


% Compressor Process - Isentropic - P as a function of v

v_Compressor = linspace(v1,v2,1000)';
T_Compressor = linspace(T1,T2,1000)';
P_Compressor = zeros(1000,1);
s_Compressor = zeros(1000,1);


for index = 1:1000
    P_Compressor(index) = C/(v_Compressor(index)^n);
    s_Compressor(index) = s1;
end

% Combustor Process - Isobaric - T as a function of s

v_Combustor = linspace(v2,v3,1000)';
T_Combustor = linspace(T2,MaxT,1000)';
P_Combustor = zeros(1000,1);
s_Combustor = zeros(1000,1);

for index = 1:1000
    s_Combustor(index) = s2 + Cp*log(T_Combustor(index)/T2);
    P_Combustor(index) = P2;
end

% Turbine Process - Isentropic - P as a function of v

v_Turbine = linspace(v3,v4,1000)';
T_Turbine = linspace(MaxT,T4,1000)';
P_Turbine = zeros(1000,1);
s_Turbine = linspace(s3,s4,1000)';

for index = 1:1000
    P_Turbine(index) = C2/(v_Turbine(index)^n);
end

% Nozzle Process - Isentropic - P as a function of v

v_Nozzle = linspace(v4,v5,1000)';
T_Nozzle = linspace(T4,T5,1000)';
P_Nozzle = zeros(1000,1);
s_Nozzle = linspace(s3,s4,1000)';

for index = 1:1000
    P_Nozzle(index) = C2/(v_Nozzle(index)^n);
end

% Reset Process - 4 to 1 - Isobaric

v_Reset = linspace(v5,v1,1000)';
P_Reset = linspace(P5,P1,1000)';
T_Reset = linspace(T5,T1,1000)';
s_Reset = zeros(1000,1);
for index = 1:1000
    P_Reset(index) = P5;
    s_Reset(index)  = s4 + Cp*log(T_Reset(index)/T5);
end

% Data Compiled for Output - 1000 points per process
Pressure_Cycle        = [P_Compressor;P_Combustor;P_Turbine;P_Nozzle;P_Reset];
SpecificVolume_Cycle  = [v_Compressor;v_Combustor;v_Turbine;v_Nozzle;v_Reset];
Temperature_Cycle     = [T_Compressor;T_Combustor;T_Turbine;T_Nozzle;T_Reset];
SpecificEntropy_Cycle = [s_Compressor;s_Combustor;s_Turbine;s_Nozzle;s_Reset];
end















