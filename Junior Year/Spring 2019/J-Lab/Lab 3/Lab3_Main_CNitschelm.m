% Lab 3
clear all
close all

%% Part 1 

R_Shunt = [55700,67700,82100,99700,149000]; % Ohms
V_Shunt = [.260,.220,.182,.150,.099]; % Volts
R_Strain = 120; % Ohms
V_Gain = V_Shunt./100; % Accounting for Amplifier

% Constant for Calculating Delta R
A = R_Strain.*R_Shunt; 
B = R_Strain + R_Shunt;

Delta_R = R_Strain - (A./B); % Equation for the change in R

[p1,s]=polyfit(Delta_R,V_Gain,1);
Fit_Data = p1(1)*Delta_R+p1(2);

Nu = length(Delta_R)-2; % First order
t_95 = tinv(0.975,Nu); % P = 95%
x_Bar = sum(Delta_R)/length(Delta_R); % Mean
Denom = sum((Delta_R-x_Bar).^2);
S_yx = (sum((V_Gain-Fit_Data).^2)/Nu).^(0.5); % Standard Error of Fit
Con_Fit = t_95*S_yx*(1/length(Delta_R)+(Delta_R-x_Bar).^2/Denom).^(0.5);
Con_Measure = t_95*S_yx*(1+1/length(Delta_R)+(Delta_R-x_Bar).^2/Denom).^(0.5);

% Plotting 
figure(1)
plot(Delta_R,V_Gain,'o',Delta_R,Fit_Data,'--',Delta_R,Fit_Data+Con_Fit,'r',Delta_R,Fit_Data-Con_Fit,'r',Delta_R,Fit_Data+Con_Measure,'b',Delta_R,Fit_Data-Con_Measure,'b');
title('Resistance vs. Amplified Output Voltage')
xlabel('\Delta R')
ylabel('\Delta e\o')
ylim([0.8*10^-3,3*10^-3])
legend('Data Points','Bestfit line','Upper Confidence Interval for Fit','Lower Confidence Interval for Fit','Upper Confidence Interval for Measurement','Lower Confidence Interval for Measurement','location','northwest') 
text(.15,1.3e-03,strcat('\Delta e\o','  ','=','  ','0.0103(Ohms)+1.851x10^-5(Volts)'))
text(.15,1.1e-03,'t-value=4.3027')
text(.15,.9e-03,'v=2')





E_i = 5; %Volts
G_F = 2.1; %Strain gauge factor
St = 2.*V_Gain./G_F*E_i;

% Strain vs. Load

Load = [.060,.110,.160,.210,.260] .* 9.81; % N
Load_Volts=[.0712-.000600,.130-.000250,.190-.000150,.248+.000200,.306+.000300]; % Volts accounting for the Zero shifting
Beam_Length = 6.58*.0254; %Length of beam in meters
Thickness = 0.050*.0254; % Thickness of beam in meters
Width = .497*.0254; % meters
E = 1.93e11; % Elastic Modulus

Strain_P=(2.*Load_Volts./100)./(G_F*E_i);
Load_Theo=(12.*Load.*(Beam_Length-Width)*(Thickness/2))./(E*(Width*Thickness^3));
[p1,s]=polyfit(Load,Strain_P,1);
Best_Fit_1=p1(1)*Load+p1(2);

figure(2)
plot(Load,Strain_P,'o',Load,Best_Fit_1,'--',Load,Load_Theo)
xlabel('Load (kg)')
ylabel('Strain (\epsilon)')
title('Applied Load vs. Strain')
legend('Data','Best Fit Line','Theoretical Data','Location','Northwest')

% Statistic Calculations 
Nu_1 = length(Load)-2; % First order
t_95_1 = tinv(0.975,Nu_1); % P = 95%
X_Bar_1 = sum(Load)/length(Load); % Mean
Denom_1 = sum((Load-X_Bar_1).^2);
s_yx_1 = (sum((Load_Theo-Best_Fit_1).^2)/Nu_1).^(0.5); % Standard Error Calculation of the Fit
Con_Fit_1 = t_95_1*s_yx_1*(1/length(Load)+(Load-X_Bar_1).^2/Denom_1).^.5;
Con_Measure_1 = t_95_1*s_yx_1*(1+1/length(Load)+(Load-X_Bar_1).^2/Denom_1).^.5;

% Plotting with Statistical Lines
figure(3)
plot(Load,Strain_P,'*',Load,Best_Fit_1,Load,Best_Fit_1+Con_Fit_1,'r',Load,Best_Fit_1-Con_Fit_1,'r',Load,Best_Fit_1+Con_Measure_1,'b',Load,Best_Fit_1-Con_Measure_1,'b')
xlabel('Force (N)')
ylabel('Strain (\epsilon)')
xlim([0.5,2.6])
title('Force vs. Strain with Confidence Intervals')
legend('Data','Bestfit','Upper Confidence Interval for Fit','Lower Confidence Interval for Fit','Upper Confidence Interval for Measurement','Lower Confidence Interval for Measurement','Location','Northwest')
text(1.5,2e-4,strcat('t','  ', '=' ,'  ', '3.184'))
text(1.5,1.5e-4,'v=3')

% Strain vs. Displacement

Disp = [0,.025,.050,.075,.100,.125,.150,.175,.200,.225,.250,.250,.225,.200,.175,.150,.125,.100,.075,.050,.025,0]*0.0254; % meters
Disp_Volts=[.0004,.0197,.0405,.0608,.0813,.1023,.1223,.1427,.1640,.1859,.2070,.2070,.1850,.1622,.1418,.1210,.0998,.0800,.0590,.0388,.0182,-.0004];
St_2=(2*Disp_Volts./100)./(G_F*E_i); % Accounting for amplifier
Disp_Theo=(3.*Disp.*(Thickness/2).*(Beam_Length-Width)./(Beam_Length^3));
[p2,s] = polyfit(Disp,St_2,1);
Best_Fit_2 = p2(1)*Disp+p2(2);




figure(4)
plot(Disp,St_2,'*',Disp,Best_Fit_2,'--',Disp,Disp_Theo)
title('Strain vs. Displacement')
xlabel('Displacement (m)')
ylabel('Strain (\epsilon)')
legend('Data','Bestfit','Theroretcial','Location','Northwest')

% Statistics
Nu_2 = length(Disp)-2; % First order
t_95_2 = tinv(0.975,Nu_2); %t value for P=95%
x_Bar_2 = sum(Disp)/length(Disp); % Mean
Denom_2 = sum((Disp-x_Bar_2).^2);
s_yx_2 = (sum((Disp_Theo-Best_Fit_2).^2)/Nu_2).^.5; %standard error of the fit
Con_Fit_2 = (t_95_2*s_yx_2*(1/length(Disp)+(Disp-x_Bar_2).^2/Denom_2).^.5);
Con_Measure_2 = (t_95_2*s_yx_2*(1+1/length(Disp)+(Disp-x_Bar_2).^2/Denom_2).^.5);

figure(5)
plot(Disp,St_2,'o',Disp,Best_Fit_2,'--',Disp,Best_Fit_2+Con_Fit_2,'r',Disp,Best_Fit_2-Con_Fit_2,'r',Disp,Best_Fit_2+Con_Measure_2,'b',Disp,Best_Fit_2-Con_Measure_2,'b')
xlabel('Displacement (m)')
ylabel('Strain (/epsilon)')
title('Displacement vs. Strain Confidence')
legend('Data','Bestfit','Upper Confidence Interval for Fit','Lower Confidence Interval for Fit','Upper Confidence Interval for Measurement','Lower Confidence Interval for Measurement','Location','Northwest')
text(4e-03,1e-4,strcat('t' , '  ','=','  ','3.1824'))
text(4e-03,.5e-4,'v=20')



% Hysteresis Calculations
Increasing = [0.0004,.0197,.0405,.0608,.0813,.1023,.1223,.1427,.1640,.1859,.2070];
Decreasing = [-.0004,.0182,.0388,.0590,.0800,.0998,.1210,.1418,.1622,.1850,.2070];
Y_Range = Increasing(end)-Increasing(1);
Difference = Increasing-Decreasing;
Hysteresis_Max = (max(abs(Difference))./(Y_Range))*100;




%% Part 2


% 2a

Header=29;
Data = importdata('Beam3.lvm','\t',Header);
Col1=Data.data(:,1);
Col2=Data.data(:,2);
for i=1:length(Col1)
    if Col1(i)>-0.005
        Base=i;
        break
    end
end

Base_Line = mean(Col2(1:Base));
Dev = std(Col2(1:Base));
Threshold = 5*Dev;

for i=1:length(Col1)
    if(abs(Col2(i)-Base_Line)>Threshold)
        Start_Time=i;
        break
    end
end

New_Time=Col1(Start_Time:length(Col1))-Col1(Start_Time);
New_Strain=Col2(Start_Time:length(Col1));
Constant=0.015;
[Locations,Values]=peakfinder(New_Strain,Constant);

figure(6)
plot(New_Time,New_Strain,New_Time(Locations),Values,'d')


% 2b
K = 3*E*Width*Thickness^3;
Denom_3=12*Beam_Length^3;
K_Beam = K/Denom_3;
Density_Steel = 7700; % kg/m^3
Mass_Hook = 0.0074; % kg
Volume_Beam=Beam_Length*Thickness*Width;
Mass_Beam = Density_Steel*Volume_Beam;
Mass_Effective = Mass_Beam/4 + Mass_Hook;
Damped_Natural_Frequency = sqrt(K_Beam/Mass_Effective); % rad/s

% 2c
Log = log(Values(2)./Values(2:length(Values)));
Con = transpose(2:length(Values));
[mb,t] = polyfit(Con,Log,1);
Best_Fit_3 = mb(1)*Con+mb(2);
alpha = mb(1);
Damp_Ratio = (alpha/sqrt(4*pi^2+alpha^2));
figure(7)
plot(Con,Log,'o',Con,Best_Fit_3)




% 2d
% Using Large Equation
for n=3:length(Values)
    Num = (1/(n-1))*log(Values(2)./Values(n));
    Damp(n) = Num./sqrt(4*pi^2+Num^2);
end

Damp_Mean = sum(Damp(3:end)./length(Damp(2:end)));
Damp_Std = std(Damp(3:length(Damp)));



%2e
k_Theo = (3*E*((Width*Thickness^3)/12))/Beam_Length^3
k_Exp = p2(1)/p1(1)

%2f
Natural_Freq_Theo = sqrt(k_Theo/Mass_Effective)
Natural_Freq_Actu = sqrt(k_Exp/Mass_Effective)



%% Part 3










































