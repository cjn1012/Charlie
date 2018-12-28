function [efficiency,Q] = Optimization(Eff_Comp,Eff_Turb,Eff_Nozz,T1,P1)


% Thermal Systems - Project 2A   
% Function to obtain performance and efficient of a jet engine with particular parameters for non constant cp

T1 = 300;
P1 = 100;
PR = 20;
MaxT = 1500;
PD=0.03;
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
v1 = (R*T1)/P1; 
s1 = interp1(Temperatures,Specific_Entropy,T1);
h1 = interp1(Specific_Entropy,Specific_Enthalpy,s1);
st1 = s1;

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
h5 = interp1(Temperatures,Specific_Enthalpy,T5);


%%%%%%
% State Calculations Non-Ideal
%%%%%%

% State 2 % Inlet of Combuster
P2n = PR*P1; 
h2n = (h2-h1)/Eff_Comp + h1;
T2n = interp1(Specific_Enthalpy,Temperatures,h2n);
st2n = interp1(Specific_Enthalpy,Specific_Entropy,h2n);
s2n  = st2n - st1 -R*log(PR) + s1;

% State 3 % Inlet of Turbine
P3n = P2; % Isobaric Process
s3n = interp1(Temperatures,Specific_Entropy,MaxT)-st2n+s2n+PD;
st3n = interp1(Temperatures,Specific_Entropy,MaxT);
h3n = interp1(Temperatures,Specific_Enthalpy,MaxT);

% State 4 % Inlet of Nozzle
h4n = h3n-(h2n-h1) ; 
h4s = ((h4n-h3n)/Eff_Turb) + h3n;
T4n = interp1(Specific_Enthalpy,Temperatures,h4n);
T4s = interp1(Specific_Enthalpy,Temperatures,h4s);
st4n = interp1(Specific_Enthalpy,Specific_Entropy,h4n);
st4s = interp1(Specific_Enthalpy,Specific_Entropy,h4s);
P4s = P3*exp((st4s-st3)/R);
P4n = P4s;
s4n = st4n - st3n - R*log(P4n/P3n) + s3n;

% State 5 % Outlet of Nozzle
P5n = P1;
s5n = s4n/Eff_Nozz;
h5n = (h5-h4n)*Eff_Nozz + h4n; 
T5n = interp1(Specific_Entropy,Temperatures,s5n);
%st5n = interp1(Specific_Enthalpy,Specific_Entropy,h5n);
%s5n = st5n - st4n - R*log(P5n/P4n) + s4n;





Q_Add    = h3n - h2n;
Q_Remove = h5n - h1;

W_Compressor = h2n - h1;
W_Turbine    = h4n - h3n;

Performance = W_Compressor + Q_Add - Q_Remove + W_Turbine;
efficiency = Performance/Q_Add;
Q = h3-h2;
end




