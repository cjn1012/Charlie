clear all; close all
clc;

Working_Fluid = {'R410a'};
Working_Fluid = Working_Fluid{1};
Max_Pressure  = 4900000;

%%%%%%%%%%%%%
% Constants for the Refridgeration Cycle for Durham, NH
%%%%%%%%%%%%%

Cold_Space = 273+20; % Kelvin - Temperature of Cold Space of Pod
Hot_Space_Liftoff  = 273 + 49; % Hottest Outside Environment
Cold_Space_Apogee   = 273 - 18; % Coldest Outside Environment

%%%%%%%%%%%%%
% Location 1 - Between Evaporator and Compressor
%%%%%%%%%%%%%

% Finding all values of the fluid at location 1
Q_1 = 1; % Saturated Vapor
T_1 = Cold_Space; % Temperature from 4 to 1
P_1 = CoolProp.PropsSI('P', 'T', T_1, 'Q', Q_1, Working_Fluid);
H_1 = CoolProp.PropsSI('H', 'T', T_1, 'Q', Q_1, Working_Fluid);
U_1 = CoolProp.PropsSI('U', 'T', T_1, 'Q', Q_1, Working_Fluid);
S_1 = CoolProp.PropsSI('S', 'T', T_1, 'Q', Q_1, Working_Fluid);
V_1 = 1/CoolProp.PropsSI('D', 'T', T_1, 'Q', Q_1, Working_Fluid);

Q_1_Apogee = 1; % Saturated Vapor
T_1_Apogee = Cold_Space; % Temperature from 4 to 1
P_1_Apogee = CoolProp.PropsSI('P', 'T', T_1_Apogee, 'Q', Q_1, Working_Fluid);
H_1_Apogee = CoolProp.PropsSI('H', 'T', T_1_Apogee, 'Q', Q_1, Working_Fluid);
U_1_Apogee = CoolProp.PropsSI('U', 'T', T_1_Apogee, 'Q', Q_1, Working_Fluid);
S_1_Apogee = CoolProp.PropsSI('S', 'T', T_1_Apogee, 'Q', Q_1, Working_Fluid);
V_1_Apogee = 1/CoolProp.PropsSI('D', 'T', T_1_Apogee, 'Q', Q_1, Working_Fluid);

%%%%%%%%%%%%%
% Location 3 - Between Condensor and Expansion Valve - Saturated Liquid
%%%%%%%%%%%%%

% Finding all values of the fluid at location 3
Q_3 = 0; % Saturated Liquid
T_3_Summer = Hot_Space_Liftoff; % Temperature from 4 to 1
P_3_Summer = CoolProp.PropsSI('P', 'T', T_3_Summer, 'Q', Q_3, Working_Fluid);
H_3_Summer = CoolProp.PropsSI('H', 'T', T_3_Summer, 'Q', Q_3, Working_Fluid);
U_3_Summer = CoolProp.PropsSI('U', 'T', T_3_Summer, 'Q', Q_3, Working_Fluid);
S_3_Summer = CoolProp.PropsSI('S', 'T', T_3_Summer, 'Q', Q_3, Working_Fluid);
V_3_Summer = 1/CoolProp.PropsSI('D', 'T', T_3_Summer, 'Q', Q_3, Working_Fluid);

T_3_Winter = Cold_Space + 5; % Temperature from 4 to 1
P_3_Winter = CoolProp.PropsSI('P', 'T', T_3_Winter, 'Q', Q_3, Working_Fluid);
H_3_Winter = CoolProp.PropsSI('H', 'T', T_3_Winter, 'Q', Q_3, Working_Fluid);
U_3_Winter = CoolProp.PropsSI('U', 'T', T_3_Winter, 'Q', Q_3, Working_Fluid);
S_3_Winter = CoolProp.PropsSI('S', 'T', T_3_Winter, 'Q', Q_3, Working_Fluid);
V_3_Winter = 1/CoolProp.PropsSI('D', 'T', T_3_Winter, 'Q', Q_3, Working_Fluid);

%%%%%%%%%%%%%
% Location 2 - Between Condensor and Expansion Valve Super-Heated Vapor
%%%%%%%%%%%%%

% Finding all values of the fluid at location 2
S_2_Summer = S_1; % Constant Entropy
P_2_Summer = P_3_Summer; % Constant Pressure Isobar
T_2_Summer = CoolProp.PropsSI('T', 'P', P_2_Summer, 'S', S_2_Summer, Working_Fluid);
H_2_Summer = CoolProp.PropsSI('H', 'P', P_2_Summer, 'S', S_2_Summer, Working_Fluid);
U_2_Summer = CoolProp.PropsSI('U', 'P', P_2_Summer, 'S', S_2_Summer, Working_Fluid);
V_2_Summer = 1/CoolProp.PropsSI('D', 'P', P_2_Summer, 'S', S_2_Summer, Working_Fluid);

P_2_Winter = P_3_Winter; % Constant Pressure Isobar
S_2_Winter = S_1_Apogee; % Constant Entropy
T_2_Winter = CoolProp.PropsSI('T', 'P', P_2_Winter, 'S', S_2_Winter, Working_Fluid);
H_2_Winter = CoolProp.PropsSI('H', 'P', P_2_Winter, 'S', S_2_Winter, Working_Fluid);
U_2_Winter = CoolProp.PropsSI('U', 'P', P_2_Winter, 'S', S_2_Winter, Working_Fluid);
V_2_Winter = 1/CoolProp.PropsSI('D', 'P', P_2_Winter, 'S', S_2_Winter, Working_Fluid);

%%%%%%%%%%%%%
% Location 2a - On Vapor Dome During the process 2 to 3 in condensor
%%%%%%%%%%%%%

% Finding all values of the fluid at location 2a
Q_2a=1; % On Vapor Dome as Saturdated Vapor
P_2a_Summer = P_2_Summer; % Isobar = Constant Pressure
T_2a_Summer = T_3_Summer; % Constant Temperature from T3
S_2a_Summer = CoolProp.PropsSI('S', 'P', P_2a_Summer, 'Q', Q_2a, Working_Fluid);
H_2a_Summer = CoolProp.PropsSI('H', 'P', P_2a_Summer, 'Q', Q_2a, Working_Fluid);
U_2a_Summer = CoolProp.PropsSI('U', 'P', P_2a_Summer, 'Q', Q_2a, Working_Fluid);
V_2a_Summer = 1/CoolProp.PropsSI('D', 'P', P_2a_Summer, 'Q', Q_2a, Working_Fluid);

P_2a_Winter = P_2_Winter; % Isobar = Constant Pressure
T_2a_Winter = T_3_Winter; % Constant Temperature from T3
S_2a_Winter = CoolProp.PropsSI('S', 'P', P_2a_Winter, 'Q', Q_2a, Working_Fluid);
H_2a_Winter = CoolProp.PropsSI('H', 'P', P_2a_Winter, 'Q', Q_2a, Working_Fluid);
U_2a_Winter = CoolProp.PropsSI('U', 'P', P_2a_Winter, 'Q', Q_2a, Working_Fluid);
V_2a_Winter = 1/CoolProp.PropsSI('D', 'P', P_2a_Winter, 'Q', Q_2a, Working_Fluid);

%%%%%%%%%%%%%
% Location 4 - Mixture between the expansion valve and evaporator
%%%%%%%%%%%%%

% Finding all values of the fluid at location 4
T_4 = T_1; % Constant Temperature from 4 to 1 Process
H_4_Summer = H_3_Summer; % Expansion Valve is ~ constant enthalpy process
P_4_Summer = P_1; % Isobaric Process
U_4_Summer = CoolProp.PropsSI('U', 'P', P_4_Summer, 'H', H_4_Summer, Working_Fluid);
S_4_Summer = CoolProp.PropsSI('S', 'P', P_4_Summer, 'H', H_4_Summer, Working_Fluid);
V_4_Summer = 1/CoolProp.PropsSI('D', 'P', P_4_Summer, 'H', H_4_Summer, Working_Fluid);

T_4_Winter = T_1_Apogee; % Constant Temperature from 4 to 1 Process
H_4_Winter = H_3_Winter; % Expansion Valve is ~ constant enthalpy process
P_4_Winter = P_1_Apogee; % Isobaric Process
U_4_Winter = CoolProp.PropsSI('U', 'P', P_4_Winter, 'H', H_4_Winter, Working_Fluid);
S_4_Winter = CoolProp.PropsSI('S', 'P', P_4_Winter, 'H', H_4_Winter, Working_Fluid);
V_4_Winter = 1/CoolProp.PropsSI('D', 'P', P_4_Winter, 'H', H_4_Winter, Working_Fluid);

%%%%%%%%%%%%%
% Calculations of the Vapor Dome for the Refridgeration Cycle
% Will be plotted as a T-s and P-h Diagram
% Two curves will be plotted for each graph, one for the Saturated Liquid and one for Saturdated Vapor Sections
%%%%%%%%%%%%%

% Constant Variables and Pressure Array to Calculate Graph Values
Q_SL = 0; % Saturated Liquid
Q_SV = 1; % Saturated Vapor
P_SL_SV = linspace(100000,Max_Pressure,1000);  % Pressures for the Saturated Liquid Curve

T_SL = zeros(length(P_SL_SV));
S_SL = zeros(length(P_SL_SV));
H_SL = zeros(length(P_SL_SV));
V_SL = zeros(length(P_SL_SV));
T_SV = zeros(length(P_SL_SV));
S_SV = zeros(length(P_SL_SV));
H_SV = zeros(length(P_SL_SV));
V_SV = zeros(length(P_SL_SV));

% Looping 1000 times to provide values for the Vapor Dome Curves for T, s, h and v. P array will be graphed with them
for index=1:1000
    T_SL(index) = CoolProp.PropsSI('T', 'P', P_SL_SV(index), 'Q', Q_SL, Working_Fluid) -273;
    S_SL(index) = CoolProp.PropsSI('S', 'P', P_SL_SV(index), 'Q', Q_SL, Working_Fluid);
    H_SL(index) = CoolProp.PropsSI('H', 'P', P_SL_SV(index), 'Q', Q_SL, Working_Fluid);
    V_SL(index) = 1/CoolProp.PropsSI('D','P',P_SL_SV(index), 'Q', Q_SL, Working_Fluid);
    T_SV(index) = CoolProp.PropsSI('T', 'P', P_SL_SV(index), 'Q', Q_SV, Working_Fluid) -273;
    S_SV(index) = CoolProp.PropsSI('S', 'P', P_SL_SV(index), 'Q', Q_SV, Working_Fluid);
    H_SV(index) = CoolProp.PropsSI('H', 'P', P_SL_SV(index), 'Q', Q_SV, Working_Fluid);
    V_SV(index) = 1/CoolProp.PropsSI('D','P',P_SL_SV(index), 'Q', Q_SV, Working_Fluid);
end

%%%%%%%%%%%%%
% Calculations of Every Point during the Refridgeration Cycle
%%%%%%%%%%%%%

% Compressor Points - Points 1 to 2 - Constant Entropy

S_Compressor_Liftoff = S_1; % Constant Entropy Process
S_Compressor_Apogee = S_1_Apogee; % Constant Entropy Process
P_Compressor_Summer = linspace(P_1,P_2_Summer,1000);
P_Compressor_Winter = linspace(P_1_Apogee,P_2_Winter,1000);
T_Compressor_Summer = zeros(length(P_Compressor_Summer));
T_Compressor_Winter = zeros(length(P_Compressor_Summer));
S_Compressor_Summer = zeros(length(P_Compressor_Summer));
S_Compressor_Winter = zeros(length(P_Compressor_Summer));
H_Compressor_Summer = zeros(length(P_Compressor_Summer));
H_Compressor_Winter = zeros(length(P_Compressor_Summer));
V_Compressor_Summer = zeros(length(P_Compressor_Summer));
V_Compressor_Winter = zeros(length(P_Compressor_Summer));

for index = 1:1000
    T_Compressor_Summer(index) = CoolProp.PropsSI('T', 'P', P_Compressor_Summer(index), 'S', S_Compressor_Liftoff, Working_Fluid) - 273;
    T_Compressor_Winter(index) = CoolProp.PropsSI('T', 'P', P_Compressor_Winter(index), 'S', S_Compressor_Apogee, Working_Fluid) - 273;
    S_Compressor_Summer(index) = CoolProp.PropsSI('S', 'P', P_Compressor_Summer(index), 'S', S_Compressor_Liftoff, Working_Fluid);
    S_Compressor_Winter(index) = CoolProp.PropsSI('S', 'P', P_Compressor_Winter(index), 'S', S_Compressor_Apogee, Working_Fluid);
    H_Compressor_Summer(index) = CoolProp.PropsSI('H', 'P', P_Compressor_Summer(index), 'S', S_Compressor_Liftoff, Working_Fluid);
    H_Compressor_Winter(index) = CoolProp.PropsSI('H', 'P', P_Compressor_Winter(index), 'S', S_Compressor_Apogee, Working_Fluid);
    V_Compressor_Summer(index) = 1/CoolProp.PropsSI('D', 'P', P_Compressor_Summer(index), 'S', S_Compressor_Liftoff, Working_Fluid);
    V_Compressor_Winter(index) = 1/CoolProp.PropsSI('D', 'P', P_Compressor_Winter(index), 'S', S_Compressor_Apogee, Working_Fluid);
end

% Condensor - Points 2 to 2a - Constant Pressure

P_Condensor_Summer = P_2_Summer; % Constant Pressure Process
P_Condensor_Winter = P_2_Winter; % Constant Pressure Process
T_Condensor_Summer = linspace(T_2_Summer,T_2a_Summer+.25,1000); % Added constant .25 so it does not go to the mixture
T_Condensor_Winter = linspace(T_2_Winter,T_2a_Winter+.25,1000); % Added constant .25 so it does not go to the mixture
S_Condensor_Summer = zeros(length(P_Condensor_Summer));
S_Condensor_Winter = zeros(length(P_Condensor_Summer));
H_Condensor_Summer = zeros(length(P_Condensor_Summer));
H_Condensor_Winter = zeros(length(P_Condensor_Summer));
V_Condensor_Summer = zeros(length(P_Condensor_Summer));
V_Condensor_Winter = zeros(length(P_Condensor_Summer));
for index = 1:1000
    S_Condensor_Summer(index) = CoolProp.PropsSI('S', 'T', T_Condensor_Summer(index), 'P', P_Condensor_Summer, Working_Fluid);
    S_Condensor_Winter(index) = CoolProp.PropsSI('S', 'T', T_Condensor_Winter(index), 'P', P_Condensor_Winter, Working_Fluid);
    H_Condensor_Summer(index) = CoolProp.PropsSI('H', 'T', T_Condensor_Summer(index), 'P', P_Condensor_Summer, Working_Fluid);
    H_Condensor_Winter(index) = CoolProp.PropsSI('H', 'T', T_Condensor_Winter(index), 'P', P_Condensor_Winter, Working_Fluid);
    V_Condensor_Summer(index) = 1/CoolProp.PropsSI('D', 'T', T_Condensor_Summer(index), 'P', P_Condensor_Summer, Working_Fluid);
    V_Condensor_Winter(index) = 1/CoolProp.PropsSI('D', 'T', T_Condensor_Winter(index), 'P', P_Condensor_Winter, Working_Fluid);
end

% Condensor - Points 2a to 3 - Constant Pressure - Straight Line inside Vapor Dome

T_Condensora_Summer = [T_2a_Summer,T_3_Summer]; % Constant Temperature
T_Condensora_Winter = [T_2a_Winter,T_3_Winter]; % Constant Temperature
P_Condensora_Summer = [P_2a_Summer,P_3_Summer];
P_Condensora_Winter = [P_2a_Winter,P_3_Winter];
S_Condensora_Summer = [S_2a_Summer,S_3_Summer];
S_Condensora_Winter = [S_2a_Winter,S_3_Winter];
H_Condensora_Summer = [H_2a_Summer,H_3_Summer];
H_Condensora_Winter = [H_2a_Winter,H_3_Winter];
V_Condensora_Summer = [V_2a_Summer,V_3_Summer];
V_Condensora_Winter = [V_2a_Winter,V_3_Winter];

% Expansion Valve - Points 3 to 4 - Constant Enthalpy

H_Valve_Summer = H_3_Summer; % Constant Enthlpy Process
H_Valve_Winter = H_3_Winter; % Constant Enthlpy Process
P_Valve_Summer = linspace(P_3_Summer,P_4_Summer,1000); % Added constant .5 so it does not go to the mixture
P_Valve_Winter = linspace(P_3_Winter,P_4_Winter,1000); % Added constant .5 so it does not go to the mixture
S_Valve_Summer = zeros(length(P_Valve_Summer));
S_Valve_Winter = zeros(length(P_Valve_Summer));
T_Valve_Summer = zeros(length(P_Valve_Summer));
T_Valve_Winter = zeros(length(P_Valve_Summer));
V_Valve_Summer = zeros(length(P_Valve_Summer));
V_Valve_Winter = zeros(length(P_Valve_Summer));
for index = 1:1000
    S_Valve_Summer(index) = CoolProp.PropsSI('S', 'H', H_Valve_Summer, 'P', P_Valve_Summer(index), Working_Fluid);
    S_Valve_Winter(index) = CoolProp.PropsSI('S', 'H', H_Valve_Winter, 'P', P_Valve_Winter(index), Working_Fluid);
    T_Valve_Summer(index) = CoolProp.PropsSI('T', 'H', H_Valve_Summer, 'P', P_Valve_Summer(index), Working_Fluid);
    T_Valve_Winter(index) = CoolProp.PropsSI('T', 'H', H_Valve_Winter, 'P', P_Valve_Winter(index), Working_Fluid);
    V_Valve_Summer(index) = 1/CoolProp.PropsSI('D', 'H', H_Valve_Summer, 'P', P_Valve_Summer(index), Working_Fluid);
    V_Valve_Winter(index) = 1/CoolProp.PropsSI('D', 'H', H_Valve_Winter, 'P', P_Valve_Winter(index), Working_Fluid);
end

% Evaporator - Points 4 to 1 - Constant Temperature - Straight Line inside Vapor Dome

T_Evaporator_Summer = [T_4,T_1]; % Constant Temperature
T_Evaporator_Winter = [T_4_Winter,T_1_Apogee]; % Constant Temperature
P_Evaporator_Summer = [P_4_Summer,P_1];
P_Evaporator_Winter = [P_4_Winter,P_1_Apogee];
S_Evaporator_Summer = [S_4_Summer,S_1];
S_Evaporator_Winter = [S_4_Winter,S_1_Apogee];
H_Evaporator_Summer = [H_4_Summer,H_1];
H_Evaporator_Winter = [H_4_Winter,H_1_Apogee];
V_Evaporator_Summer = [V_4_Summer,V_1];
V_Evaporator_Winter = [V_4_Winter,V_1_Apogee];

%%%%%%%%%%%%%
% Graphing the T-s and P-h Graphs with Vapor Dome and Labels
%%%%%%%%%%%%%

% T-s Graph for the Summer
figure(1)
% Vapor Dome
plot(S_SL/1000,T_SL,'k',S_SV/1000,T_SV,'k')
hold on
% Processes

% Winter
plot(S_Compressor_Winter/1000,T_Compressor_Winter,'b',S_Condensor_Winter/1000, T_Condensor_Winter-273.15,'b',S_Condensora_Winter/1000, T_Condensora_Winter-273.15,'b',S_Valve_Winter/1000,T_Valve_Winter-273.15,'b',S_Evaporator_Winter/1000, T_Evaporator_Winter-273.15,'b')
text(S_1_Apogee/1000,T_1_Apogee-273.15, '\leftarrow State 1')
text(S_2_Winter/1000,T_2_Winter-273.15, '\leftarrow State 2')
text(S_2a_Winter/1000-.025,T_2a_Winter-273.15-1, 'State 2a','FontSize', 8)
text(S_3_Winter/1000-.04,T_3_Winter-273.15, 'State 3')
text(S_4_Winter/1000-.02,T_4_Winter-273.15-1, 'State 4')

% Summer
plot(S_Compressor_Summer/1000,T_Compressor_Summer,'r',S_Condensor_Summer/1000, T_Condensor_Summer-273.15,'r',S_Condensora_Summer/1000, T_Condensora_Summer-273.15,'r',S_Valve_Summer/1000,T_Valve_Summer-273.15,'r',S_Evaporator_Summer/1000, T_Evaporator_Summer-273.15,'m')
plot(S_Compressor_Winter/1000,T_Compressor_Winter,'m')
text(S_2_Summer/1000,T_2_Summer-273.15, '\leftarrow State 2')
text(S_2a_Summer/1000,T_2a_Summer-273.15, '\leftarrow State 2a')
text(S_3_Summer/1000-.04,T_3_Summer-273.15, 'State 3')
text(S_4_Summer/1000-.01,T_4-273.15-1, 'State 4')

title('Refrigeration Cycle during the Summer in Durham, NH')
xlabel('Entropy (KJ/K)')
ylabel('Temperature (Celcius)')
xlim([.8 2])
ylim([-20 80])
hold off

% Syntax
%title('Refridgeration Cycle during the Winter in Durham, NH','FontSize',20)
xlabel('Entropy (J/K)','FontSize',22)
set(gca,'fontsize',20)
ylabel('Temperature (Celcius)','FontSize',22)
set(gca,'fontsize',20)
lgd = legend('\color{red} Summer','\color{blue} Winter','\color{black} Vapor Dome');
lgd.FontSize = 22;
xlim([.9 1.9])
ylim([0 75])
hold off

% P-h Graph for the Summer
figure(2)
% Vapor Dome
plot(H_SL/1000,P_SL_SV/1000,'k',H_SV/1000,P_SL_SV/1000,'k')
hold on

% Processes

%Winter
P_Condensor_Winter = [P_2_Winter,P_2a_Winter];
H_Condensor_Winter = [H_2_Winter,H_2a_Winter];
P_Valve_Winter = [P_3_Winter,P_4_Winter];
H_Valve_Winter = [H_3_Winter,H_4_Winter];
plot(H_Compressor_Winter/1000,P_Compressor_Winter/1000,'b',H_Condensor_Winter/1000, P_Condensor_Winter/1000,'b',H_Condensora_Winter/1000, P_Condensora_Winter/1000,'b',H_Valve_Winter/1000,P_Valve_Winter/1000,'b',H_Evaporator_Winter/1000, P_Evaporator_Winter/1000,'b')
text(H_1_Apogee/1000,P_1_Apogee/1000, '\leftarrow State 1')
text(H_2_Winter/1000,P_2_Winter/1000, '\leftarrow State 2')
text(H_3_Winter/1000-12,P_3_Winter/1000, 'State 3')
text(H_4_Winter/1000-5,P_4_Winter/1000-35, 'State 4')

% Summer
P_Condensor_Summer = [P_2_Summer,P_2a_Summer];
H_Condensor_Summer = [H_2_Summer,H_2a_Summer];
P_Valve_Summer = [P_3_Summer,P_4_Summer];
H_Valve_Summer = [H_3_Summer,H_4_Summer];
plot(H_Compressor_Summer/1000,P_Compressor_Summer/1000,'r',H_Condensor_Summer/1000, P_Condensor_Summer/1000,'r',H_Condensora_Summer/1000, P_Condensora_Summer/1000,'r',H_Valve_Summer/1000,P_Valve_Summer/1000,'r',H_Evaporator_Summer/1000, P_Evaporator_Summer/1000,'m')
plot(H_Compressor_Winter/1000,P_Compressor_Winter/1000,'m')
text(H_2_Summer/1000,P_2_Summer/1000, '\leftarrow State 2')
text(H_3_Summer/1000-12,P_3_Summer/1000, 'State 3')
text(H_4_Summer/1000-5,P_4_Summer/1000-50, 'State 4')

% Plot Syntax
%title('Refrigeration Cycle during the Summer and Winter in Durham, NH','FontSize',20)
xlabel('Enthalpy (KJ/K)','FontSize',22)
set(gca,'fontsize',20)
ylabel('Pressure (KPa)','FontSize',22)
set(gca,'fontsize',20)
lgd = legend('\color{red} Summer','\color{blue} Winter','\color{black} Vapor Dome');
lgd.FontSize = 22;
xlim([150 475])
ylim([500 5000])
hold off

% P-h Graph for the Summer
figure(3)
% Vapor Dome
plot(V_SL/1000,P_SL_SV/1000,'k',V_SV/1000,P_SL_SV/1000,'k')
hold on

% Processes

%Winter
P_Condensor_Winter = [P_2_Winter,P_2a_Winter];
V_Condensor_Winter = [V_2_Winter,V_2a_Winter];
P_Valve_Winter = [P_3_Winter,P_4_Winter];
V_Valve_Winter = [V_3_Winter,V_4_Winter];
plot(V_Compressor_Winter/1000,P_Compressor_Winter/1000,'b',V_Condensor_Winter/1000, P_Condensor_Winter/1000,'b',V_Condensora_Winter/1000, P_Condensora_Winter/1000,'b',V_Valve_Winter/1000,P_Valve_Winter/1000,'b',V_Evaporator_Winter/1000, P_Evaporator_Winter/1000,'b')
text(V_1/1000,P_1/1000+20, '\leftarrow State 1')
text(V_2_Winter/1000,P_2_Winter/1000+20, '\leftarrow State 2')
text(V_3_Winter/1000+.00000035,P_3_Winter/1000-35, 'State 3')
text(V_4_Winter/1000,P_4_Winter/1000-35, 'State 4')

% Summer
P_Condensor_Summer = [P_2_Summer,P_2a_Summer];
V_Condensor_Summer = [V_2_Summer,V_2a_Summer];
P_Valve_Summer = [P_3_Summer,P_4_Summer];
V_Valve_Summer = [V_3_Summer,V_4_Summer];
plot(V_Compressor_Summer/1000,P_Compressor_Summer/1000,'r',V_Condensor_Summer/1000, P_Condensor_Summer/1000,'r',V_Condensora_Summer/1000, P_Condensora_Summer/1000,'r',V_Valve_Summer/1000,P_Valve_Summer/1000,'r',V_Evaporator_Summer/1000, P_Evaporator_Summer/1000,'m')
plot(V_Compressor_Winter/1000,P_Compressor_Winter/1000,'m')
text(V_2_Summer/1000,P_2_Summer/1000+20, '\leftarrow State 2')
text(V_3_Summer/1000+.0000005,P_3_Summer/1000-65, 'State 3')
text(V_4_Summer/1000,P_4_Summer/1000-70, 'State 4')
text(275,2500, 'Cycle Starts at State 1 towards State 2','FontSize',16)

% Plot Syntax
%title('Refrigeration Cycle during the Summer and Winter in Durham, NH','FontSize',20)
xlabel('Specific Volume (m3/Kg)','FontSize',22)
set(gca,'fontsize',20)
ylabel('Pressure (KPa)','FontSize',22)
set(gca,'fontsize',20)
lgd = legend('\color{red} Summer','\color{blue} Winter','\color{black} Vapor Dome');
lgd.FontSize = 22;
xlim([.0000005 .000025])
ylim([1000 5000])
hold off

%%%%%%%%%%%%%
% COP vs Outside Temperature for 3 Different Refrigerants
%%%%%%%%%%%%%

Working_Fluid = {'R410a'};
Working_Fluid_2 = {'R290'};
Working_Fluid_3 = {'Ammonia'};
Working_Fluids = [Working_Fluid,Working_Fluid_2,Working_Fluid_3];

Outside_Temperature = linspace(273 - 18, 273 + 49, 1000);

COP_R410a   = zeros(1000);
COP_R290    = zeros(1000);
COP_Ammonia = zeros(1000);

j=1;
for i = Working_Fluids
    k=1;
    for temps = Outside_Temperature
        
        Working_Fluid = Working_Fluids{j};
        %%%%%%%%%%%%%
        % Constants for the Refridgeration Cycle for Durham, NH
        %%%%%%%%%%%%%
        
        
        Cold_Space = 273+20; % Kelvin - Temperature of Cold Space of Pod
        
        %Hot_Space_Liftoff  = 273 + 49; % Hottest Outside Environment
        %Cold_Space_Apogee   = 273 - 18; % Coldest Outside Environment
        
        
        %%%%%%%%%%%%%
        % Location 1 - Between Evaporator and Compressor
        %%%%%%%%%%%%%
        
        % Finding all values of the fluid at location 1
        Q_1 = 1; % Saturated Vapor
        T_1 = Cold_Space; % Temperature from 4 to 1
        P_1 = CoolProp.PropsSI('P', 'T', T_1, 'Q', Q_1, Working_Fluid);
        H_1 = CoolProp.PropsSI('H', 'T', T_1, 'Q', Q_1, Working_Fluid);
        U_1 = CoolProp.PropsSI('U', 'T', T_1, 'Q', Q_1, Working_Fluid);
        S_1 = CoolProp.PropsSI('S', 'T', T_1, 'Q', Q_1, Working_Fluid);
        V_1 = 1/CoolProp.PropsSI('D', 'T', T_1, 'Q', Q_1, Working_Fluid);
        
        
        
        %%%%%%%%%%%%%
        % Location 3 - Between Condensor and Expansion Valve - Saturated Liquid
        %%%%%%%%%%%%%
        
        % Finding all values of the fluid at location 3
        Q_3 = 0; % Saturated Liquid
        T_3_Summer = temps + 5; % Temperature from 4 to 1
        P_3_Summer = CoolProp.PropsSI('P', 'T', T_3_Summer, 'Q', Q_3, Working_Fluid);
        H_3_Summer = CoolProp.PropsSI('H', 'T', T_3_Summer, 'Q', Q_3, Working_Fluid);
        U_3_Summer = CoolProp.PropsSI('U', 'T', T_3_Summer, 'Q', Q_3, Working_Fluid);
        S_3_Summer = CoolProp.PropsSI('S', 'T', T_3_Summer, 'Q', Q_3, Working_Fluid);
        V_3_Summer = 1/CoolProp.PropsSI('D', 'T', T_3_Summer, 'Q', Q_3, Working_Fluid);
        
        
        
        
        %%%%%%%%%%%%%
        % Location 2 - Between Condensor and Expansion Valve Super-Heated Vapor
        %%%%%%%%%%%%%
        
        % Finding all values of the fluid at location 2
        S_2_Summer = S_1; % Constant Entropy
        P_2_Summer = P_3_Summer; % Constant Pressure Isobar
        T_2_Summer = CoolProp.PropsSI('T', 'P', P_2_Summer, 'S', S_2_Summer, Working_Fluid);
        H_2_Summer = CoolProp.PropsSI('H', 'P', P_2_Summer, 'S', S_2_Summer, Working_Fluid);
        U_2_Summer = CoolProp.PropsSI('U', 'P', P_2_Summer, 'S', S_2_Summer, Working_Fluid);
        V_2_Summer = 1/CoolProp.PropsSI('D', 'P', P_2_Summer, 'S', S_2_Summer, Working_Fluid);
        
        
        %%%%%%%%%%%%%
        % Location 2a - On Vapor Dome During the process 2 to 3 in condensor
        %%%%%%%%%%%%%
        
        % Finding all values of the fluid at location 2a
        Q_2a=1; % On Vapor Dome as Saturdated Vapor
        P_2a_Summer = P_2_Summer; % Isobar = Constant Pressure
        T_2a_Summer = T_3_Summer; % Constant Temperature from T3
        S_2a_Summer = CoolProp.PropsSI('S', 'P', P_2a_Summer, 'Q', Q_2a, Working_Fluid);
        H_2a_Summer = CoolProp.PropsSI('H', 'P', P_2a_Summer, 'Q', Q_2a, Working_Fluid);
        U_2a_Summer = CoolProp.PropsSI('U', 'P', P_2a_Summer, 'Q', Q_2a, Working_Fluid);
        V_2a_Summer = 1/CoolProp.PropsSI('D', 'P', P_2a_Summer, 'Q', Q_2a, Working_Fluid);
        
        
        %%%%%%%%%%%%%
        % Location 4 - Mixture between the expansion valve and evaporator
        %%%%%%%%%%%%%
        
        % Finding all values of the fluid at location 4
        T_4 = T_1; % Constant Temperature from 4 to 1 Process
        H_4_Summer = H_3_Summer; % Expansion Valve is ~ constant enthalpy process
        P_4_Summer = P_1; % Isobaric Process
        U_4_Summer = CoolProp.PropsSI('U', 'P', P_4_Summer, 'H', H_4_Summer, Working_Fluid);
        S_4_Summer = CoolProp.PropsSI('S', 'P', P_4_Summer, 'H', H_4_Summer, Working_Fluid);
        V_4_Summer = 1/CoolProp.PropsSI('D', 'P', P_4_Summer, 'H', H_4_Summer, Working_Fluid);
        
        
        
        
        % Co-Efficient of Performance where COP = (QL/Wnet) = ((h1-h4)/(h2-h1))
        
        COP_Summer = (H_1 - H_4_Summer)/(H_2_Summer - H_1);
        
        % Mass Flow Rate where flow rate = Cooling Capacity / Q_L
        
        dH_Evaporator_Summer = H_1 - H_4_Summer;
        Cooling_Capacity = 5000; % Watts
        Mass_Flow_Rate_Summer = Cooling_Capacity / dH_Evaporator_Summer;
        
        % Compressor Power
        
        Compressor_Power_Summer = Mass_Flow_Rate_Summer * (H_2_Summer - H_1); % [J/s]
        Final_Data_Names  = ['COP_Summer',' Mass_Flow_Rate_Summer',' Compressor_Power_Summer'];
        Final_Data_Values = [COP_Summer,Mass_Flow_Rate_Summer,Compressor_Power_Summer];
        
        if j == 1
            COP_R410a(k)   = COP_Summer;
        elseif j == 2
            COP_R290(k)   = COP_Summer;
        else
            COP_Ammonia(k)   = COP_Summer;
        end
        
        k=k+1;
    end
    j=j+1;
end

for pos = 1:1000
    COP_R410a(pos) = abs(COP_R410a(pos));
    COP_R290(pos) = abs(COP_R290(pos));
    COP_Ammonia(pos) = abs(COP_Ammonia(pos));
end

% Plotting COP vs Outside Temperature
figure(4)
semilogx(COP_R410a,Outside_Temperature-273,'r',COP_R290,Outside_Temperature-273,'b',COP_Ammonia,Outside_Temperature-273,'k')

% Plot Syntax
%title('COP of a Thermal Mangement System with Varying Outside Temperature','FontSize',18)
xlabel('Coefficient of Performance','FontSize',22)
set(gca,'fontsize',20)
ylabel('Temperature','FontSize',22)
set(gca,'fontsize',20)
xlim([5 1000])
ylim([-20 50])
lgd = legend('\color{red} R410a','\color{blue} R290','\color{black} Ammonia');
lgd.FontSize = 22;
hold off

