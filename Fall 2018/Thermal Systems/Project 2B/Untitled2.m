function Performance = Performance(PR,MaxT)


% Thermal Systems - Project 2A   
% Function to obtain performance and efficient of a jet engine with particular parameters for non constant cp

T1 = 300;
P1 = 100;
Eff_Comp = .72;
Eff_Turb = .98;
Eff_Nozz = .90;

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
st2n = st1 + R*log(PR);
s2n  = interp1(Specific_Enthalpy,Specific_Entropy,h2n)-st2+s1;

% State 3 % Inlet of Turbine
P3n = P2-300; % Isobaric Process
s3n = interp1(Temperatures,Specific_Entropy,MaxT)-st2n+s2n;
st3n = interp1(Temperatures,Specific_Entropy,MaxT);
h3n = interp1(Temperatures,Specific_Enthalpy,MaxT);

% State 4 % Inlet of Nozzle
h4n = -Eff_Turb * (h3n-h4) + h3n ; 
st4n = interp1(Specific_Enthalpy,Specific_Entropy,h4n);
s4n = interp1(Specific_Enthalpy,Specific_Entropy,h4n)-st3n+s3n;

T4n = interp1(Specific_Enthalpy,Temperatures,h4n);
P4n = P4;


% State 5 % Outlet of Nozzle
h5n = -Eff_Nozz * (h4n-h5) + h4n ; 
s5n = interp1(Specific_Enthalpy,Specific_Entropy,h5n);
st5n = interp1(Specific_Enthalpy,Specific_Entropy,h5n);
T5n = interp1(Specific_Enthalpy,Temperatures,h5n);
P5n = P1;


Q_Add    = h3n - h2n;
Q_Remove = h5n - h1;

W_Compressor = h2n - h1;
W_Turbine    = h4n - h3n;

Performance = W_Compressor + Q_Add - Q_Remove + W_Turbine;
efficiency = Performance/Q_Add;

end




