clear all; close all
clc;


%%%%%%%%%%%%%
% COP vs Outside Temperature for 3 Different Refrigerants
%%%%%%%%%%%%%

Working_Fluid_1  = {'R410a'};
Working_Fluid_2  = {'R290'};
Working_Fluid_3  = {'Ammonia'};
Working_Fluid_4 = {'Acetone'};

Working_Fluids   = [Working_Fluid_1,Working_Fluid_2,Working_Fluid_3,Working_Fluid_4];

Outside_Temperature = linspace(273 - 22, 273 + 52, 1000);

COP_R410a       = zeros(1000,1);
COP_R290        = zeros(1000,1);
COP_Ammonia     = zeros(1000,1);
COP_Acetone     = zeros(1000,1);

MFR_R410a       = zeros(1000,1);
MFR_R290        = zeros(1000,1);
MFR_Ammonia     = zeros(1000,1);
MFR_Acetone     = zeros(1000,1);


MFRW_R410a       = zeros(1000,1);
MFRW_R290        = zeros(1000,1);
MFRW_Ammonia     = zeros(1000,1);
MFRW_Acetone     = zeros(1000,1);

Comp_Power_R410a =zeros(1000,1);
Comp_Power_R290 =zeros(1000,1);
Comp_Power_Ammonia = zeros(1000,1);
Comp_Power_Acetone =zeros(1000,1);

Comp_PowerW_R410a =zeros(1000,1);
Comp_PowerW_R290 =zeros(1000,1);
Comp_PowerW_Ammonia = zeros(1000,1);
Comp_PowerW_Acetone =zeros(1000,1);
j=1;
for i = Working_Fluids
    k=1;
    for temps = Outside_Temperature
        
        Working_Fluid_1 = Working_Fluids{j};
        %%%%%%%%%%%%%
        % Constants for the Refridgeration Cycle for Durham, NH
        %%%%%%%%%%%%%
        
        Cold_Space = 273+20; % Kelvin - Temperature of Cold Space of Pod
        
        %%%%%%%%%%%%%
        % Location 1 - Between Evaporator and Compressor
        %%%%%%%%%%%%%
        
        % Finding all values of the fluid at location 1
        Q_1 = 1; % Saturated Vapor
        T_1 = Cold_Space; % Temperature from 4 to 1
        P_1 = CoolProp.PropsSI('P', 'T', T_1, 'Q', Q_1, Working_Fluid_1);
        H_1 = CoolProp.PropsSI('H', 'T', T_1, 'Q', Q_1, Working_Fluid_1);
        U_1 = CoolProp.PropsSI('U', 'T', T_1, 'Q', Q_1, Working_Fluid_1);
        S_1 = CoolProp.PropsSI('S', 'T', T_1, 'Q', Q_1, Working_Fluid_1);
        V_1 = 1/CoolProp.PropsSI('D', 'T', T_1, 'Q', Q_1, Working_Fluid_1);
        
        %%%%%%%%%%%%%
        % Location 3 - Between Condensor and Expansion Valve - Saturated Liquid
        %%%%%%%%%%%%%
        
        % Finding all values of the fluid at location 3
        Q_3 = 0; % Saturated Liquid
        T_3_Summer = temps + 5; % Temperature from 4 to 1
        P_3_Summer = CoolProp.PropsSI('P', 'T', T_3_Summer, 'Q', Q_3, Working_Fluid_1);
        H_3_Summer = CoolProp.PropsSI('H', 'T', T_3_Summer, 'Q', Q_3, Working_Fluid_1);
        U_3_Summer = CoolProp.PropsSI('U', 'T', T_3_Summer, 'Q', Q_3, Working_Fluid_1);
        S_3_Summer = CoolProp.PropsSI('S', 'T', T_3_Summer, 'Q', Q_3, Working_Fluid_1);
        V_3_Summer = 1/CoolProp.PropsSI('D', 'T', T_3_Summer, 'Q', Q_3, Working_Fluid_1);
        

        %%%%%%%%%%%%%
        % Location 2 - Between Condensor and Expansion Valve Super-Heated Vapor
        %%%%%%%%%%%%%
        
        % Finding all values of the fluid at location 2
        S_2_Summer = S_1; % Constant Entropy
        P_2_Summer = P_3_Summer; % Constant Pressure Isobar
        T_2_Summer = CoolProp.PropsSI('T', 'P', P_2_Summer, 'S', S_2_Summer, Working_Fluid_1);
        H_2_Summer = CoolProp.PropsSI('H', 'P', P_2_Summer, 'S', S_2_Summer, Working_Fluid_1);
        U_2_Summer = CoolProp.PropsSI('U', 'P', P_2_Summer, 'S', S_2_Summer, Working_Fluid_1);
        V_2_Summer = 1/CoolProp.PropsSI('D', 'P', P_2_Summer, 'S', S_2_Summer, Working_Fluid_1);
        
        %%%%%%%%%%%%%
        % Location 2a - On Vapor Dome During the process 2 to 3 in condensor
        %%%%%%%%%%%%%
        
        % Finding all values of the fluid at location 2a
        Q_2a=1; % On Vapor Dome as Saturdated Vapor
        P_2a_Summer = P_2_Summer; % Isobar = Constant Pressure
        T_2a_Summer = T_3_Summer; % Constant Temperature from T3
        S_2a_Summer = CoolProp.PropsSI('S', 'P', P_2a_Summer, 'Q', Q_2a, Working_Fluid_1);
        H_2a_Summer = CoolProp.PropsSI('H', 'P', P_2a_Summer, 'Q', Q_2a, Working_Fluid_1);
        U_2a_Summer = CoolProp.PropsSI('U', 'P', P_2a_Summer, 'Q', Q_2a, Working_Fluid_1);
        V_2a_Summer = 1/CoolProp.PropsSI('D', 'P', P_2a_Summer, 'Q', Q_2a, Working_Fluid_1);
        
        %%%%%%%%%%%%%
        % Location 4 - Mixture between the expansion valve and evaporator
        %%%%%%%%%%%%%
        
        % Finding all values of the fluid at location 4
        T_4 = T_1; % Constant Temperature from 4 to 1 Process
        H_4_Summer = H_3_Summer; % Expansion Valve is ~ constant enthalpy process
        P_4_Summer = P_1; % Isobaric Process
        U_4_Summer = CoolProp.PropsSI('U', 'P', P_4_Summer, 'H', H_4_Summer, Working_Fluid_1);
        S_4_Summer = CoolProp.PropsSI('S', 'P', P_4_Summer, 'H', H_4_Summer, Working_Fluid_1);
        V_4_Summer = 1/CoolProp.PropsSI('D', 'P', P_4_Summer, 'H', H_4_Summer, Working_Fluid_1);
        
        % Co-Efficient of Performance where COP = (QL/Wnet) = ((h1-h4)/(h2-h1))
        
        COP_Summer = (H_1 - H_4_Summer)/(H_2_Summer - H_1);
        
        % Mass Flow Rate where flow rate = Cooling Capacity / Q_L
        
        dH_Evaporator_Summer = H_1 - H_4_Summer;
        Cooling_Capacity = 500+14000; % Watts
        Cooling_Capacity_Winter = 2000-13500; % Watts
        Mass_Flow_Rate_Summer = Cooling_Capacity / dH_Evaporator_Summer;
        Mass_Flow_Rate_Winter = Cooling_Capacity_Winter / dH_Evaporator_Summer;
        
        % Compressor Power
        
        Compressor_Power_Summer = Mass_Flow_Rate_Summer * (H_2_Summer - H_1); % [J/s]
        Compressor_Power_Winter = Mass_Flow_Rate_Winter * (H_2_Summer - H_1); % [J/s]
      
        if j == 1
            COP_R410a(k)   = COP_Summer;
            MFR_R410a(k)   = Mass_Flow_Rate_Summer;
            MFRW_R410a(k)  = Mass_Flow_Rate_Winter;
            Comp_Power_R410a(k) = Compressor_Power_Summer;
            Comp_PowerW_R410a(k) = Compressor_Power_Winter;
        elseif j == 2
            COP_R290(k)   = COP_Summer;
            MFR_R290(k)   = Mass_Flow_Rate_Summer;
            MFRW_R290(k)   = Mass_Flow_Rate_Winter;
            Comp_Power_R290(k) = Compressor_Power_Summer;
            Comp_PowerW_R290(k) = Compressor_Power_Winter;
        elseif j == 3
            COP_Ammonia(k)   = COP_Summer;
            MFR_Ammonia(k)   = Mass_Flow_Rate_Summer;
            MFRW_Ammonia(k)   = Mass_Flow_Rate_Winter;
            Comp_Power_Ammonia(k) = Compressor_Power_Summer;
            Comp_PowerW_Ammonia(k) = Compressor_Power_Winter;
        else
            COP_Acetone(k)   = COP_Summer;
            MFR_Acetone(k)   = Mass_Flow_Rate_Summer;
            MFRW_Acetone(k)   = Mass_Flow_Rate_Winter;
            Comp_Power_Acetone(k) = Compressor_Power_Summer;
            Comp_PowerW_Acetone(k) = Compressor_Power_Winter;
        end
        k=k+1;
    end
    j=j+1;
end

for pos = 1:1000
    COP_R410a(pos) = abs(COP_R410a(pos));
    COP_R290(pos) = abs(COP_R290(pos));
    COP_Ammonia(pos) = abs(COP_Ammonia(pos));
    COP_Acetone(pos) = abs(COP_Acetone(pos));
end

% Create a list from best to worst refrigerants in respect to its COP, MFR, Compressor Power and Q dot

COP_R410a_Max       = COP_R410a(973); % The index selection for Takeoff Temperature
COP_R290_Max        = COP_R290(973);
COP_Ammonia_Max     = COP_Ammonia(973);
COP_Acetone_Max     = COP_Acetone(973);

COPW_R410a_Max       = COP_R410a(28); % The index selection for Takeoff Temperature
COPW_R290_Max        = COP_R290(28);
COPW_Ammonia_Max     = COP_Ammonia(28);
COPW_Acetone_Max     = COP_Acetone(28);

Comp_R410a_Max       = Comp_Power_R410a(973); % The index selection for Takeoff Temperature
Comp_R290_Max        = Comp_Power_R290(973);
Comp_Ammonia_Max     = Comp_Power_Ammonia(973);
Comp_Acetone_Max     = Comp_Power_Acetone(973);

CompW_R410a_Max       = Comp_PowerW_R410a(28); % The index selection for Takeoff Temperature
CompW_R290_Max        = Comp_PowerW_R290(28);
CompW_Ammonia_Max     = Comp_PowerW_Ammonia(28);
CompW_Acetone_Max     = Comp_PowerW_Acetone(28);

MFR_R410a_Max       = MFR_R410a(973); % The index selection for Takeoff Temperature
MFR_R290_Max        = MFR_R290(973);
MFR_Ammonia_Max     = MFR_Ammonia(973);
MFR_Acetone_Max     = MFR_Acetone(973);

MFRW_R410a_Max       = MFRW_R410a(28); % The index selection for Takeoff Temperature
MFRW_R290_Max        = MFRW_R290(28);
MFRW_Ammonia_Max     = MFRW_Ammonia(28);
MFRW_Acetone_Max     = MFRW_Acetone(28);

ACOP_TakeOff_Sorted =sortrows({COP_R410a_Max,'R410a';COP_R290_Max,'R290';COP_Ammonia_Max,'Ammonia';COP_Acetone_Max,'Acetone'},1);
ACOPW_TakeOff_Sorted =sortrows({COPW_R410a_Max,'R410a';COPW_R290_Max,'R290';COPW_Ammonia_Max,'Ammonia';COPW_Acetone_Max,'Acetone'},1);
AComp_TakeOff_Sorted =sortrows({Comp_R410a_Max,'R410a';Comp_R290_Max,'R290';Comp_Ammonia_Max,'Ammonia';Comp_Acetone_Max,'Acetone'},1);
ACompW_TakeOff_Sorted =sortrows({CompW_R410a_Max,'R410a';CompW_R290_Max,'R290';CompW_Ammonia_Max,'Ammonia';CompW_Acetone_Max,'Acetone'},1);
AMFR_TakeOff_Sorted =sortrows({MFR_R410a_Max,'R410a';MFR_R290_Max,'R290';MFR_Ammonia_Max,'Ammonia';MFR_Acetone_Max,'Acetone'},1);
AMFRW_TakeOff_Sorted =sortrows({MFRW_R410a_Max,'R410a';MFRW_R290_Max,'R290';MFRW_Ammonia_Max,'Ammonia';MFRW_Acetone_Max,'Acetone'},1);

% Plotting COP vs Outside Temperature
figure(1)
% Plot Lines
hold on
semilogx(COP_R410a,Outside_Temperature-273,'r',COP_R290,Outside_Temperature-273,'b')
semilogx(COP_Acetone,Outside_Temperature-273,'g',COP_Ammonia,Outside_Temperature-273,'m')
plot([5,100],[-20,-20],[5,150],[50,50])
hold off
% Plot Syntax
text(9,48,'Takeoff Environment')
text(12,-18,'Apogee Environment')
xlabel('Coefficient of Performance','FontSize',22)
set(gca,'fontsize',20)
ylabel('Temperature','FontSize',22)
set(gca,'fontsize',20)
xlim([5 100])
ylim([-25 55])
lgd = legend('\color{red} R410a','\color{blue} R290','\color{green} Acetone','\color{magenta} Ammonia');
lgd.FontSize = 14;
hold off


% Plot of Outside Temperature vs Mass Flow Rate with Lines for extreme environment temperatures
figure(2)
hold on
plot(Outside_Temperature-273,MFR_R410a,'r',Outside_Temperature-273,MFR_R290,'b',Outside_Temperature-273,MFR_Ammonia,'g',Outside_Temperature-273,MFR_Acetone,'m')
plot(Outside_Temperature-273,MFRW_R410a,'r',Outside_Temperature-273,MFRW_R290,'b',Outside_Temperature-273,MFRW_Ammonia,'g',Outside_Temperature-273,MFRW_Acetone,'m')
plot([50,50],[0,.6])
plot([-20,-20],[-.4,0])
plot([-25,60],[0,0],'k')
% Plot Syntax
text(0,.03,'Cooling to Heating Line','Fontsize',8)
text(12,.3,'Takeoff Environment')
text(-19,-.3,'Apogee Environment')
xlabel('Temperature (C)','FontSize',22)
set(gca,'fontsize',20)
ylabel('Mass Flow Rate (kg/s)','FontSize',22)
ylim([-.4,.6])
xlim([-25,60])
set(gca,'fontsize',20)
lgd = legend('\color{red} R410a','\color{blue} R290','\color{green} Acetone','\color{magenta} Ammonia','Location','northwest');
lgd.FontSize = 10;
hold off


% Plot of Outside Temperature vs required compressor power with extreme temperature lines graphed
figure(3)
hold on
plot(Outside_Temperature-273,Comp_Power_R410a,'r',Outside_Temperature-273,Comp_Power_R290,'b',Outside_Temperature-273,Comp_Power_Ammonia,'g',Outside_Temperature-273,Comp_Power_Acetone,'m')
plot(Outside_Temperature-273,Comp_PowerW_R410a,'r',Outside_Temperature-273,Comp_PowerW_R290,'b',Outside_Temperature-273,Comp_PowerW_Ammonia,'g',Outside_Temperature-273,Comp_PowerW_Acetone,'m')
plot([50,50],[0,15000])
plot([-20,-20],[-10000,0])
plot([-25,60],[0,0],'k')
% Plot Syntax
text(-20,500,'Cooling to Heating Line','Fontsize',8)
text(13,10000,'Takeoff Environment')
text(-19,-1600,'Apogee Environment')
xlabel('Temperature (C)','FontSize',22)
set(gca,'fontsize',20)
ylabel('Compressor Power (J/s)','FontSize',22)
ylim([-10000,15000])
xlim([-25,60])
set(gca,'fontsize',20)
lgd = legend('\color{red} R410a','\color{blue} R290','\color{green} Acetone','\color{magenta} Ammonia','Location','northwest');
lgd.FontSize = 10;
hold off
