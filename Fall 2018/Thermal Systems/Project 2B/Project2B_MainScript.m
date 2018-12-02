%% Part 1
% eff = .6016
% A study of the effects of compressor, turbine and nozzle efficiencity to the overall cycle

% Sum of efficiencies cannot exceed 260%
    % Nozzle can not exceed 98%
    % Turbine can not exceed 90%
    % Compressor can not exceed 90%
    
% Choose the best percentages and graph a T-s diagram

clear all
close all

CompressorEff = linspace(.70,.90,21);
TurbineEff = linspace(.78,.98,21);
NozzleEff = linspace(.70,.90,21);

EffMatrix = zeros(200,4);
l = 1;

for i = CompressorEff
    for j = TurbineEff
        for k = NozzleEff
            if i+j+k == 2.60
                EffMatrix(l,1) = i;
                EffMatrix(l,2) = j;
                EffMatrix(l,3) = k;
                EffMatrix(l,4) = Optimization(i,j,k);
                l = l + 1;
            end
        end
    end
end

                
Sorted = sort(EffMatrix,4);                

Eff_Comp = Sorted(1,1);
Eff_Turb = Sorted(1,2);
Eff_Nozz = Sorted(1,3);

[EntropyIdeal,TemperatureIdeal,EntropyActual,TemperatureActual,h3,h2] = T_s_Data(Eff_Comp,Eff_Turb,Eff_Nozz);

figure(1)
plot(EntropyIdeal,TemperatureIdeal,'b')
hold on
plot(EntropyActual,TemperatureActual,'r')

Q_Comb = h3-h2 % = 930

%% Part 2

% Study the effects of pressure ratio and maximum temperature on the actual cycle.
% ? Show plots representing:
% – the efficiency as a function of the pressure ratio and maximum temperature
% – the increase in kinetic energy as a function of the pressure ratio and maximum
% temperature


PR = linspace(1,50,50);
MaxT = linspace(1500,2500,50);
Efficiency_PR = zeros(50,1);
Efficiency_MaxT = zeros(50,1);
Performance_PR = zeros(50,1);
Performance_MaxT = zeros(50,1);
x=1
for index = PR
    [Performance_PR(x),Efficiency_PR(x)] = Efficiency(index,1500);
    x=x+1
end
x=1
for index = MaxT
    [Performance_MaxT(x),Efficiency_MaxT(x)] = Efficiency(30,index);
    x=x+1
end

figure(2)
plot(PR,Efficiency_PR)

figure(3)
plot(MaxT,Efficiency_MaxT)

figure(4)
plot(PR,Performance_PR)

Figure(5)
plot(MaxT,Performance_MaxT)




%% Part 3

% Estimate the range of the plane as a function of its speed.
% ? Assume the aircraft behaves as a giant NACA 0015 airfoil [in terms of lift and drag
% characteristics], and select a reasonable plan-form area of the aircraft. Assume that
% Qadded?cycle ? Wdrag. The mass of the UAV is expected to be approximately 1500
% kg [without fuel and equipment]. Assuming a fixed amount of fuel that can be stored
% aboard, evaluate the ranges of the aircraft as a function of cruising speed.
Vv = linspace(100,1000,100);
range2 = zeros(100,1)';
x=1
for index = Vv
    range2(x) = Range(index) ;
    x = x+1
end


plot(Vv,range2)
    
    





%% Part 4

% Recommend a fuel and evaluate CO2 emissions.
% ? Study at least three different fuels.
% ? Assuming a fixed amount of fuel that can be stored aboard, evaluate the ranges of the
% aircraft as a function of cruising speed for each fuel.
% ? Compare the CO2 emissions of the drone for the 3 fuels.

% The gases - Kerosene, Kerosene gasoline, Avgas
















































