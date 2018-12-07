%% Part 1
% eff = .6016
% A study of the effects of compressor, turbine and nozzle efficiencity to the overall cycle

% Sum of efficiencies cannot exceed 260%
    % Nozzle can not exceed 98%
    % Turbine can not exceed 90%
    % Compressor can not exceed 90%
    
% Choose the best percentages and graph a T-s diagram

CompressorEff = linspace(.70,.90,21);
TurbineEff = linspace(.70,.90,21);
NozzleEff = linspace(.78,.98,21);
EffMatrix = zeros(200,4);
l = 1;

for i = CompressorEff
    for j = TurbineEff
        for k = NozzleEff
            if i+j+k == 2.60
                EffMatrix(l,1) = i;
                EffMatrix(l,2) = j; 
                EffMatrix(l,3) = k;
                EffMatrix(l,4) = Optimization(i,j,k,300,100);
                l = l + 1;
            end
        end
    end
end

Sorted = sortrows(EffMatrix,4);                

Eff_Comp = Sorted(end,1);
Eff_Turb = Sorted(end,2);
Eff_Nozz = Sorted(end,3);
[EntropyIdeal,TemperatureIdeal,EntropyActual,TemperatureActual,h3,h2] = T_s_Data(Eff_Comp,Eff_Turb,Eff_Nozz);
Q = h3-h2;

figure(1)
plot(EntropyIdeal,TemperatureIdeal,'b',EntropyActual,TemperatureActual,'r')
hold on
plot(EntropyIdeal(3001:4000),TemperatureIdeal(3001:4000),'m')

text(EntropyActual(1)+.01   ,TemperatureActual(1)   , '\leftarrow Compressor Inlet','FontSize',12)
text(EntropyActual(1001)+.01,TemperatureActual(1001)-7, '\leftarrow Combustor Inlet','FontSize',12)
text(EntropyActual(2001)+.01,TemperatureActual(2001), '\leftarrow Turbine Inlet','FontSize',12)
text(EntropyActual(3001)+.01,TemperatureActual(3001), '\leftarrow Nozzle Inlet','FontSize',12)
text(EntropyActual(4001)+.01,TemperatureActual(4001), '\leftarrow Nozzle Outlet','FontSize',12)
text(1.55,1100, 'Q = 885 KJ/kg','FontSize',12)
annotation('textarrow',[.26,.28],[.25,.4],'String','Cycle Direction','FontSize',10)

xlabel('Specific Entropy [KJ/kg]','FontSize',20)
set(gca,'fontsize',18)
ylabel('Temperature [K]','FontSize',20)
set(gca,'fontsize',18)
lgd = legend('\color{blue} Ideal Cycle','Location','northwest','\color{red} Actual Cycle');
lgd.FontSize = 14;
xlim([1.25,3.5])
ylim([0,1600])
hold off


%% Part 2

% Study the effects of pressure ratio and maximum temperature on the actual cycle.
% ? Show plots representing:
% – the efficiency as a function of the pressure ratio and maximum temperature
% – the increase in kinetic energy as a function of the pressure ratio and maximum
% temperature

PR = linspace(1,50,50);
MaxT = linspace(1000,2000,50);
Efficiency_PR = zeros(50,1);
Efficiency_MaxT = zeros(50,1);
Performance_PR = zeros(50,1);
Performance_MaxT = zeros(50,1);
x=1;
for index = PR
    [Performance_PR(x),Efficiency_PR(x)] = Efficiency(index,1500,100,300);
    x=x+1;
end
x=1;
for index = MaxT
    [Performance_MaxT(x),Efficiency_MaxT(x)] = Efficiency(20,index,100,300);
    x=x+1;
end


figure(2)

yyaxis left
plot(PR,Efficiency_PR,'b')
ylabel('Efficiency [ ]','FontSize',20,'Color','b')
set(gca,'fontsize',18,'fontsize',18)
ylim([0,.7])

yyaxis right
plot(PR,Performance_PR,'r')
ylabel('Performance [KJ/kg]','FontSize',20,'Color','r')
set(gca,'fontsize',18,'fontsize',18)
ylim([0,600])

lgd = legend('\color{blue} Efficiency','\color{red} Performance','Location','northwest');
lgd.FontSize = 14;
xlabel('Pressure Ratio [KJ/kg]','FontSize',20)
xlim([0,50])
hold off


figure(3)

yyaxis left
plot(MaxT,Efficiency_MaxT,'b')
ylabel('Efficiency [ ]','FontSize',20,'Color','b')
set(gca,'fontsize',18,'fontsize',18)
ylim([.48,.58])

yyaxis right
plot(MaxT,Performance_MaxT,'r')
ylabel('Performance [KJ/kg]','FontSize',20,'Color','r')
set(gca,'fontsize',18,'fontsize',18)
ylim([0,800])

lgd = legend('\color{blue} Efficiency','\color{red} Performance','Location','northwest');
lgd.FontSize = 14;
xlabel('Maximum Temperature [K]','FontSize',20)
xlim([1000,2000])
hold off

%% Part 3

% Estimate the range of the plane as a function of its speed.
% ? Assume the aircraft behaves as a giant NACA 0015 airfoil [in terms of lift and drag
% characteristics], and select a reasonable plan-form area of the aircraft. Assume that
% Qadded?cycle ? Wdrag. The mass of the UAV is expected to be approximately 1500
% kg [without fuel and equipment]. Assuming a fixed amount of fuel that can be stored
% aboard, evaluate the ranges of the aircraft as a function of cruising speed.

Vv = linspace(0,200,100);
range_kero = zeros(100,1)';
Q_Kero = 42000000;

x=1;
for index = Vv
    range_kero(x) = Range(index,.9,Q_Kero) ;
    x = x+1;
end


figure(4)
plot(Vv,range_kero)
xlabel('Velocity [m/s]','FontSize',20)
set(gca,'fontsize',18)
ylabel('Range [km]','FontSize',20)
set(gca,'fontsize',18)
xlim([30,200])
ylim([0,10000])

text(35,400, '\leftarrow 30 m/s are required for flight (Weight = Lift)','FontSize',12)
text(41,8850, '\leftarrow 39 m/s yields the best range overall','FontSize',12)
hold off

%% Part 4

% Recommend a fuel and evaluate CO2 emissions.
% ? Study at least three different fuels.
% ? Assuming a fixed amount of fuel that can be stored aboard, evaluate the ranges of the
% aircraft as a function of cruising speed for each fuel.
% ? Compare the CO2 emissions of the drone for the 3 fuels.

% The gases - Kerosene, Kerosene gasoline, Avgas


Vv = linspace(30,200,100);
range_kero = zeros(100,1)';
range_methanol = zeros(100,1)';
range_TNT = zeros(100,1)';
Q_kero = 42000000;
Q_methanol = 20000000;
Q_TNT = 4600000;
Q_Fuels = [Q_kero,Q_methanol,Q_TNT];

x=1;
for index = Vv
    range_kero(x) = Range(index,.9,Q_Fuels(1)) ;
    x = x+1;
end

x=1;
for index = Vv
    range_methanol(x) = Range(index,.9,Q_Fuels(2)) ;
    x = x+1;
end

x=1;
for index = Vv
    range_TNT(x) = Range(index,.9,Q_Fuels(3)) ;
    x = x+1;
end


figure(5)
plot(Vv,range_kero,'r',Vv,range_methanol,'m',Vv,range_TNT,'b')
xlabel('Velocity [m/s]','FontSize',20)
set(gca,'fontsize',18)
ylabel('Range [km]','FontSize',20)
set(gca,'fontsize',18)
lgd = legend('\color{red} Kerosene','\color{magenta} Methanol','\color{blue} TNT','Location','northeast');
lgd.FontSize = 14;
xlim([30,200])
ylim([0,10000])

text(100,6000, 'Q(Kerosene) = 42,000,000 J/kg','FontSize',12)
text(100,5000, 'Q(Methanol) = 20,000,000 J/kg','FontSize',12)
text(100,4000, 'Q(TNT)      = 4,600,000 J/kg','FontSize',12)

hold off




%% Part 5


Air_Data = xlsread('air_data1.xls');

Q_kero = 42000000;
Q_methanol = 19700000;
Q_TNT = 4610000;
Q_Fuels = [Q_kero,Q_methanol,Q_TNT];

Altitudes = Air_Data(1:43,6);
Densities = Air_Data(1:43,10);
Pressures = Air_Data(1:43,8);
Temperatures = Air_Data(1:43,7);
range_altitudes = zeros(43,1)';
performance_altitudes = zeros(43,1);
for index = 1:43
    range_altitudes(index) = Range(75,Densities(index),Q_Fuels(1));
    performance_altitudes(index) = Efficiency(30,1500,Pressures(index),Temperatures(index)+273);
end

figure(6)
yyaxis left
plot(Altitudes(1:43),range_altitudes(1:43),'b')
ylabel('Range [km]','FontSize',20,'Color','b')
set(gca,'fontsize',18,'fontsize',18)
ylim([2000,10000])


yyaxis right
plot(Altitudes(1:43),performance_altitudes(1:43),'r')
ylabel('Performance [KJ/kg]','FontSize',20,'Color','r')
set(gca,'fontsize',18,'fontsize',18)
ylim([450,600])

lgd = legend('\color{blue} Range','\color{red} Performance','Location','northwest');
lgd.FontSize = 14;
xlim([0,18000])
text(14200,579,'Best Case','fontsize',12)

xlabel('Altitude [m]','FontSize',20)


%% Part 6
% Additional task. Perform an additional tasks of similar complexity to tasks 1, 3 or 4.
% ? The task should require analysis using Matlab.



Altitudes = Air_Data(1:43,6);
Densities = Air_Data(1:43,10);
Masses = linspace(50,500,43);
Range2 = zeros(5,5);
Performance = zeros(5,5);
x=1
y=1

for index = 1:9:43
    for index2 = 1:9:43
        [Range2(x,y),Performance(x,y) ]= RangeExtra(100,Densities(index),43000000,Masses(index2));
        y=y+1
    end
    x=x+1
end


















