
function Pressure_Cycle,SpecificVolume_Cycle = CycleDataPv(PR,MaxT)

% Thermal Systems - Project 2A   
% Task 1

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
s1 = 0;
C=(P1*(v1^n));  % Isentropic Constant
 
% State 2 % Inlet of Combuster  
P2 = PR*P1; 
v2 = (C/P2)^(1/n);   
T2= (P2*v2)/R;  
s2 = 0;
% State 3 % Inlet of Turbine   
P3 = P2; % Isobaric Process
v3 = (R * MaxT)/P3; 
s3 = Cp*log(MaxT/T2);
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

v_Compressor = linspace(v1,v2,1000);
P_Compressor = zeros(1000,1);
s_Compressor = zeros(1000,1);
T_Compressor = zeros(1000,1);

for index = 1:1000
    P_Compressor(index) = C/(v_Compressor(index)^n);
end

% Combustor Process - Isobaric - T as a function of s

T_Combustor = linspace(T2,MaxT,1000)';
s_Combustor = zeros(1000,1);
P_Combustor = zeros(1000,1);
v_Combustor = linspace(v2,v3,1000);

for index = 1:1000
    s_Combustor(index) = Cp*log(T_Combustor(index)/T2);
    P_Combustor(index) = P2;
end

% Turbine Process - Isentropic - P as a function of v

v_Turbine = linspace(v3,v4,1000);
P_Turbine = zeros(1000,1);

for index = 1:1000
    P_Turbine(index) = C2/(v_Turbine(index)^n);
end


% Nozzle Process - Isentropic - P as a function of v

v_Nozzle = linspace(v4,v5,1000);
P_Nozzle = zeros(1000,1);

for index = 1:1000
    P_Nozzle(index) = C2/(v_Nozzle(index)^n);
end


% Reset Process - 4 to 1 - Isobaric

v_Reset = linspace(v5,v1,1000);
P_Reset = linspace(P5,P1,1000)';


Pressure_Cycle = [P_Compressor;P_Combustor;P_Turbine;P_Nozzle;P_Reset];
SpecificVolume_Cycle = [v_Compressor,v_Combustor,v_Turbine,v_Nozzle,v_Reset];

plot(SpecificVolume_Cycle,Pressure_Cycle)

%plot(v_Compressor,P_Compressor,v_Combustor,P_Combustor,v_Turbine,P_Turbine,v_Nozzle,P_Nozzle,v_Reset,P_Reset)
%hold on
%plot(s_Combustor,T_Combustor)













