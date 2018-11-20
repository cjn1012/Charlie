% A study of the effects of compressor, turbine and nozzle efficiencity to the overall cycle

% Sum of efficiencies cannot exceed 260%
    % Nozzle can not exceed 98%
    % Turbine can not exceed 90%
    % Compressor can not exceed 90%
    
% Choose the best percentages and graph a T-s diagram

%%%%%%
% Definition of Constants
%%%%%%

Air_Data = xlsread('air_data1.xls');
Temperatures    = Air_Data(:,1);
Specific_Enthalpy = Air_Data(:,3);
Specific_Entropy = Air_Data(:,4);

Efficiency_Compressor = .9;
Efficiency_Turbine    = .95;
Efficiency_Nozzle     = .90;

R = 0.287;     % Constant
PR = 10;
MaxT =1500;

% Ideal Case

% State 1 % Inlet of Compressor 
T1a = 300;
P1a = 100;
s1a = interp1(Temperatures,Specific_Entropy,T1a);
h1a = interp1(Specific_Entropy,Specific_Enthalpy,s1a);
st1a = s1a;

% State 2 % Inlet of Combuster  
P2a = PR*P1a; 
s2a = s1a;
st2a = st1a + R*log(PR);
T2a = interp1(Specific_Entropy,Temperatures,s2a);
h2a = interp1(Temperatures,Specific_Enthalpy,T2a);

% State 3 % Inlet of Turbine
P3a = P2a; % Isobaric Process
s3a = interp1(Temperatures,Specific_Entropy,MaxT)-st2a+s2a;
st3a = interp1(Temperatures,Specific_Entropy,MaxT);
h3a = interp1(Temperatures,Specific_Enthalpy,MaxT);

% State 4 % Inlet of Nozzle
s4a = s3;
T4a = interp1(Specific_Entropy,Temperatures,s4a);
st4a = interp1(Temperatures,Specific_Entropy,T4a);
h4a = h3a - (h2a-h1a);
P4a = P3a*(exp((st4a-st3a)/R));

% State 5 % Outlet of Nozzle
s5a = s4a;
P5a = P1a; 
st5a = st4a + R*log(P5a/P4a);
T5a = interp1(Specific_Entropy,Temperatures,st5a) ;
h5a = interp1(Temperatures,Specific_Enthalpy,T5a);

%%%%%%
% State Calculations for Actual Cycle
%%%%%%

% State 1 % Inlet of Compressor 
s1 = s1a;
h1 = h1a;
st1 = st1a;

% State 2 % Inlet of Combuster  
P2 = P2a;
st2 = st1 + R*log(PR);
h2 = ((h2a-h1)/Efficiency_Compressor)+h1;
T2 = interp1(Specific_Enthalpy,Temperatures,h2);
s2 = interp1(Specific_Enthalpy,Specific_Entropy,h2)-st1+s1;

% State 3 % Inlet of Turbine
P3 = P2; % Isobaric Process
s3 = interp1(Temperatures,Specific_Entropy,MaxT)-st2+s2;
st3 = interp1(Temperatures,Specific_Entropy,MaxT);
h3 = interp1(Temperatures,Specific_Enthalpy,MaxT);

% State 4 % Inlet of Nozzle
h4 = h3 - Efficiency_Turbine*(h3-h4a);
T4 = interp1(Specific_Enthalpy,Temperatures,h4);
s4 = interp1(Temperatures,Specific_Entropy,T4)-st3+s3;
st4 = interp1(Temperatures,Specific_Entropy,T4);

% State 5 % Outlet of Nozzle
P5 = P1; 
h5 = h4 - Efficiency_Nozzle*(h4-h5a);
s5 = interp1(Specific_Enthalpy,Specific_Entropy,h5);
T5 = interp1(Specific_Enthalpy,Temperatures,h5);


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


% Compressor Process - Isentropic - P as a function of v


s_Compressor = linspace(s1,s2,1000)';
T_Compressor = linspace(T1,T2,1000)';


% Combustor Process - Isobaric - T as a function of s

T_Combustor = linspace(T2,MaxT,1000)';
s_Combustor = zeros(1000,1);

for index = 1:1000
    s_Combustor(index) = (interp1(Temperatures,Specific_Entropy,T_Combustor(index)))-st2+s2;
end

% Turbine Process - Isentropic - P as a function of v

s_Turbine = linspace(s3,s4,1000)';
T_Turbine = linspace(MaxT,T4,1000)';



% Nozzle Process - Isentropic - P as a function of v

s_Nozzle = linspace(s4,s5,1000)';
T_Nozzle = linspace(T4,T5,1000)';

% Reset Process - 4 to 1 - Isobaric

T_Reset = linspace(T5,T1,1000)';
s_Reset = zeros(1000,1);

for index = 1:1000
    s_Reset(index) = interp1(Temperatures,Specific_Entropy,T_Reset(index))-st1+s1;
end

% Data Compiled for Output - 1000 points per process
Temperature_Cycle     = [T_Compressor;T_Combustor;T_Turbine;T_Nozzle;T_Reset];
SpecificEntropy_Cycle = [s_Compressor;s_Combustor;s_Turbine;s_Nozzle;s_Reset];

figure(1)
plot(s_Compressor,T_Compressor)
figure(2)
plot(SpecificEntropy_Cycle,Temperature_Cycle)
























