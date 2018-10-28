clear all;
close all;

% Worst case scenario. Delta T constant
% Consider the effect of temperature difference between the cycle and the environments on the length, diameter, and number of
% pipes, the flow speeds in the tubes, and the mass of the heat exchanger. Consider at least two
% different pipe materials. 

% Provide plots and optionally tables to show all dependencies. Specify
% the design parameters of each of the two heat exchangers, and the expected temperature
% differences between the cycle and the environments. 

% Worst Case Scenario is at takeoff, refridgerator
Temperature_Environement = 273 + 50; % Environement Temperature at Takeoff
Temperature_Pod = 273 + 20; % Environement Temperature at Takeoff
Length = 1; % Meters

Q_Takeoff_Fiberglass = 7300; 

Convection_Coeff_R290 = (120)/2000;
Convection_Coeff_Air = 17;
Convection_Coeff_Copper = 398;
Convection_Coeff_PVC = .136;

Q_Copper = zeros(1000,1);
Q_PVC = zeros(1000,1);
i=1;
r = linspace(.01,5,1000);
for r1 = r
        
        Q_Copper(i) = (2*pi*(Temperature_Environement-Temperature_Pod)*(Length))/((1/(Convection_Coeff_R290*r1/.0008))+(log(1+(.0008/r1))/Convection_Coeff_Copper)+(1/(Convection_Coeff_Air/.0008*(r1+.0008))));
        Q_PVC(i) = (2*pi*(Temperature_Environement-Temperature_Pod)*(Length))/((1/(Convection_Coeff_R290*r1/.0008))+(log(1+(.0008/r1))/Convection_Coeff_PVC)+(1/(Convection_Coeff_Air/.0008*(r1+.0008))));
        i = i+1;
end

j=1;
for i = Q_Copper'
    if (i > Q_Takeoff_Fiberglass) == 1
        x=j;
    else
        
        j= j+1;
    end
end
for i = Q_PVC'
    if (i > Q_Takeoff_Fiberglass) == 1
        y=j;
    else
        j= j+1;
    end
end




mdot = .0455 %kg/s
R_total_Copper = r(x);
R_total_PVC = r(y);
Num_Pipes_Copper = R_total_Copper/.03/1.5 % Length of pipe equals 1.5
v_Copper = mdot/(8960*pi*.03^2);
volume_Copper = pi*((.03+.0008)^2-.03^2)*1.5*Num_Pipes_Copper;
mass_Copper = 8960*volume_Copper

Num_Pipes_PVC = R_total_PVC/.03/1.5 % Length of pipe equals 1.5
v_PVC = mdot/(1380*pi*.03^2);
volume_PVC = pi*((.03+.0008)^2-.03^2)*1.5*Num_Pipes_PVC;
mass_PVC = 1380*volume_PVC

Num_Pipes = r/.03/1.5; % Length of pipe equals 1.5;

figure(1)
plot(Q_Copper,Num_Pipes,'r')
hold on
plot(Q_PVC,Num_Pipes,'b')
xlim([0,14000])
ylim([0,30])
plot([7300,7300],[0,30],'k')
set(gca,'fontsize',18)
xlabel('Heat Transfer for Cooling','FontSize',20)
ylabel('Number of 1.5 Meter Pipes','FontSize',20)
lgd = legend('\color{red} Copper Piping','\color{blue} PVC Piping','Location','northwest');
lgd.FontSize = 10;
text(7400,5,'Required Heat Transfer for Takeoff','Color','k','FontSize',9)
text(7400,3.5,'With a Mass Flow Rate = .0455 kg/s','Color','k','FontSize',9)
text(7400,11,strcat('Mass of Copper = ', num2str(mass_Copper)),'Color','red','FontSize',10)
text(3000,17,strcat('Mass of PVC = ', num2str(mass_PVC)),'Color','blue','FontSize',10)


