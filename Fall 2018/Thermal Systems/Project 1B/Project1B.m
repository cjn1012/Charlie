clear all; close all
clc;


Apogee_Temperature  = 260; % Kelvin
Takeoff_Temperature = 300; % Kelvin

% The array of best and worst case scenarios
Temperatures = [Apogee_Temperature,Takeoff_Temperature];
T_Pod = 288; % Kelvin
% Convection coefficient for difference in temperatures
h_Convection_Coefficient_Air = [35,17]; %W/m2k

% Pod Dimensions
Length_Pod = 2; % Meters
Width_Pod  = 1; % Meters
Height_Pod = 1; % Meters

% Material Densities
Aluminum_Density        = 2700; %kg/m3, 6061
Fiberglass_Density      = 2550; %kg/m3
Magnesium_AZ61_Density = 1800; %kg/m3
Magnesium_AZ31_Density = 1800; %kg/m3
Densities = [Aluminum_Density,Fiberglass_Density,Magnesium_AZ61_Density,Magnesium_AZ31_Density]; % 4 Total, loop every two

% Thermal Conductivity for 260 and 300 Kelvin
Aluminum_Thermal_k       = [147,155];      %W/mk
Fiberglass_Thermal_k     = [0.0325,0.038]; %W/mk
Magnesium_AZ61_Thermal_k = [82,85];        %W/mk
Magnesium_AZ31_Thermal_k = [60,63.5];      %W/mk
Thermal_k_Array_T1 = [Aluminum_Thermal_k(1),Fiberglass_Thermal_k(1),Magnesium_AZ61_Thermal_k(1),Magnesium_AZ31_Thermal_k(1)]; % 4 Total
Thermal_k_Array_T2 = [Aluminum_Thermal_k(2),Fiberglass_Thermal_k(2),Magnesium_AZ61_Thermal_k(2),Magnesium_AZ31_Thermal_k(2)]; % 4 Total


% Wall thickness array
% Wall_Thickness = linspace(.001,.01,10);

Total_Mass_Pod = zeros(4,10);
R_Pod_Total    = zeros(4,10);
Q_Dot_Environment = zeros(4,10);
i=1;
for density = Densities
    j=1;
    for Wall_Thickness = linspace(.001,.01,10)
        Temperature = Temperatures(1);
        h_Convection_Coefficient_Air_T1 = h_Convection_Coefficient_Air(1);
        Pod_Surface_Area_1 = 2*Height_Pod*Width_Pod;
        Pod_Surface_Area_2 = 3*Length_Pod*Height_Pod; % Not including Top Face against Drone
        Total_Mass_Pod(i,j) = Wall_Thickness*density * (Pod_Surface_Area_1 + Pod_Surface_Area_2); % Total Weight of the Pod Frame
        % Resistances
        R_Conductivity_Pod_1 = Wall_Thickness / (Thermal_k_Array_T1(i) * Pod_Surface_Area_1);
        R_Conductivity_Pod_2 = Wall_Thickness / (Thermal_k_Array_T1(i) * Pod_Surface_Area_2);
        R_Convection_Air_1 = 1/(h_Convection_Coefficient_Air_T1*Pod_Surface_Area_1);
        R_Convection_Air_2 = 1/(h_Convection_Coefficient_Air_T1*Pod_Surface_Area_2);
        R_Pod_Total(i,j) = 1/(1/(R_Conductivity_Pod_1 + R_Convection_Air_1) + 1/(R_Conductivity_Pod_2+ R_Convection_Air_2));
        Q_Dot_Environment(i,j) = T_Pod-Temperature/(R_Pod_Total(i,j));
        j = j+1;
    end
    i = i+1;
end
    
Total_Mass_Pod_T2 = zeros(4,10);
R_Pod_Total_T2    = zeros(4,10);
Q_Dot_Environment_T2 = zeros(4,10);
i=1;
for density = Densities
    j=1;
    for Wall_Thickness = linspace(.001,.01,10)
        
        Temperature = Temperatures(2);
        h_Convection_Coefficient_Air_T2 = h_Convection_Coefficient_Air(2);
        Pod_Surface_Area_1 = 2*Height_Pod*Width_Pod;
        Pod_Surface_Area_2 = 3*Length_Pod*Height_Pod; % Not including Top Face against Drone
        Total_Mass_Pod_T2(i,j) = Wall_Thickness.*density * (Pod_Surface_Area_1 + Pod_Surface_Area_2); % Total Weight of the Pod Frame
        % Resistances
        R_Conductivity_Pod_1 = Wall_Thickness ./ (Thermal_k_Array_T2(i) * Pod_Surface_Area_1);
        R_Conductivity_Pod_2 = Wall_Thickness ./ (Thermal_k_Array_T2(i) * Pod_Surface_Area_2);
        R_Convection_Air_1 = 1/(h_Convection_Coefficient_Air_T2*Pod_Surface_Area_1);
        R_Convection_Air_2 = 1/(h_Convection_Coefficient_Air_T2*Pod_Surface_Area_2);
        R_Pod_Total_T2(i,j) = 1/(1/(R_Conductivity_Pod_1 + R_Convection_Air_1) + 1/(R_Conductivity_Pod_2+ R_Convection_Air_2));
    
        Q_Dot_Environment_T2(i,j) = T_Pod-Temperature/(R_Pod_Total_T2(i,j));
        j = j+1;
    end
    i = i+1;
end



plot(Total_Mass_Pod')
ylabel('Kilograms')
xlabel('thickness mm')

plot(Q_Dot_Environment([1,3,4],1:10)')


