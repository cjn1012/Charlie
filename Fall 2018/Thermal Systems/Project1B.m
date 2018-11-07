  
% Thermal Systems - Project 2A   
% Task 1

%%%%%%
% Definition of Constants
%%%%%%
R=0.287;     % Constant
n=1.4;       % Constant
k = (1/0.4); % Constant
r=28;        %Pressure Ratio 
Cp = 1.005; % kJ/Kg

%%%%%%
% State Calculations
%%%%%%

% State 1 % Inlet of Compressor 
T1 = 300;       % Inlet Temperature
P1=100;         % Inlet Pressure
v1 = (R*T1)/P1; 

C=(P1*(v1^n));  % Isentropic Constant
 
% State 2 % Inlet of Combuster  
P2 = r*P1; 
v2 = (C/P2)^(1/n);   
T2= (P2*v2)/R;  
     
% State 3 % Inlet of Turbine   
P3 = P2; % Isobaric Process
T3 = 1540 + 273;  
v3 = (R * T3)/P3; 
   
% State 4 % Inlet of Nozzle
T4 = T3 - (T2-T1);   
v4 = ((P3 * (v3^n))/(R*T4))^(k); 
P4 = (R * T4)/v4;

% State 5 % Outlet of Nozzle
P5 = P1; 
v5 = ((P4*(v4^n))/P5)^(1/n); 
T5 = (P5*v5)/R;   


%%%%%%
% Process Calculations
%%%%%%

% Compressor - Isentropic
% Combustor  - Isobaric
% Turbine    - Isentropic
% Nozzle     - Isentropic


% Compressor Process - Isentropic - P as a function of v

v_Compressor = linspace(v1,v2,1000);
P_Compressor = zeros(1000,1);

for index = 1:1000
    P_Compressor(index) = C/(v_Compressor(index)^n);
end

% Combustor Process - Isobaric - T as a function of s

T_Combustor = linspace(T2,T3,1000);
s_Combustor = zeros(1000,1);



















