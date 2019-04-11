clear all 
close all

% Part 1
Temp_I = 21.2+273.15; % Deg C
Temp_F = 22.6+273.15; % Deg C
Pres_I = 102268; % Pa
Pres_F = 102235; % Pa
R = 8.314; % J/molK
MW = 28.97/1000; % g/mol

Density_I = (MW*Pres_I)/(R*Temp_I); % kg/m3
Density_F = (MW*Pres_F)/(R*Temp_F); % kg/m3

% Part 2
AoA_Foil_0  = [1.6,3.2,2.95,3.54,3.85,3.57,3.79,3.65,3.52,3.05,3.31,3.16,2.98,2.55,3.65,3.44,3.21,3.5,3.4,3.31,2.91,3.14,.49].*0.0254;
AoA_Foil_9  = [1.6,3.2,7,5.9,5.25,4.4,4.2,3.8,3.55,3.05,3.25,3.15,1.6,2.2,2.6,2.65,2.7,3.05,3.1,3.05,2.7,2.9,.45].*0.0254;
AoA_Foil_18 = [1.6,3.25,4.05,3.95,4.00,3.7,4.1,4.15,4.2,3.65,4.0,3.85,1.5,1.9,2.25,2.35,2.5,2.9,3.0,3.1,2.85,3.25,0.5].*0.0254;

Den_Water = 998; %kg/m3
Gravity = 9.81; %m/s2

Pres_Foil_0 =  Pres_I - Den_Water*Gravity.*AoA_Foil_0;
Pres_Foil_9 =  Pres_I - Den_Water*Gravity.*AoA_Foil_9;
Pres_Foil_18 =  Pres_I - Den_Water*Gravity.*AoA_Foil_18;

% Part 3
Vel_Foil_0 = sqrt((2/Density_I)*(Pres_Foil_0(1)-Pres_Foil_0(2)));
Vel_Foil_9 = sqrt((2/Density_I)*(Pres_Foil_9(1)-Pres_Foil_9(2)));
Vel_Foil_18 = sqrt((2/Density_I)*(Pres_Foil_18(1)-Pres_Foil_18(2)));
Vel_Foil_Mean = mean([Vel_Foil_0,Vel_Foil_9,Vel_Foil_18]);

% Part 4
AoAs = [-9,-6,-3,0,3,6,9,12,15,18,21];
L_Force = [-9.3,-6.3,-3,0.4,3.3,6.3,9,11.5,13.5,11.1,10.8];
D_Force = [1.5,1.1,0.8,0.75,0.85,1.1,1.5,1.95,2.5,5.1,5.7];
A = .15*.30; % Chord * Span
C_L = L_Force/(.5*Density_I*(Vel_Foil_Mean^2)*A);
C_D = D_Force/(.5*Density_I*(Vel_Foil_Mean^2)*A);

C_L_Fit = polyfit(AoAs,C_L,1);
Zero = -C_L_Fit(2)/C_L_Fit(1); % This says we dont need to adjust at all, below +-0.5 degrees



% Part 5

plot(AoAs,C_L,AoAs,C_D)
xlabel('Angle of Attack (\circ)')
ylabel('Coefficient of Lift/Drag ( )')
legend('Coefficent of Lift','Coefficient of Drag','location','northwest')


% Part 6

plot(AoAs,C_L./C_D)
xlabel('Angle of Attack (\circ)')
ylabel('Ratio of C_L and C_D ( )')


% Part 7

AoA_Foil_0_1  = AoA_Foil_0(3:22)-(AoA_Foil_0(end)*.0254);
AoA_Foil_9_1  = AoA_Foil_9(3:22)-(AoA_Foil_9(end)*.0254);
AoA_Foil_18_1 = AoA_Foil_18(3:22)-(AoA_Foil_18(end)*.0254);

Pres_Foil_0_1 =  -Den_Water*Gravity.*AoA_Foil_0_1;
Pres_Foil_9_1 =  -Den_Water*Gravity.*AoA_Foil_9_1;
Pres_Foil_18_1 = -Den_Water*Gravity.*AoA_Foil_18_1;

c = .15; % Chorde Length
t = .12;% Thickness to chorde
L = .30; % Span

x_top = [0.76, 3.81, 11.43, 19.05, 38, 62, 80.77, 101.35, 121.92, 137.16]./1000;
x_bot = [1.52,7.62,15.24,22.86,41.15,59.44,77.73,96.02,114.3,129.54]./1000; 

y_top =5*t.*( (.2969.*sqrt(x_top/c))-(.1260.*x_top/c)-(.3516.*(x_top/c).^2)+(.2843.*(x_top/c).^3)-(.1015.*(x_top/c).^4));
y_bot =-5*t.*( (.2969.*sqrt(x_bot/c))-(.1260.*x_bot/c)-(.3516.*(x_bot/c).^2)+(.2843.*(x_bot/c).^3)-(.1015.*(x_bot/c).^4));

plot(x_top,y_top,x_bot,y_bot)


for i = 2:length(y_top)
    theta = atan((y_top(i)-y_top(i-1))/(x_top(i)-x_top(i-1)));
    F_Top_x(i-1) = L*(x_top(i)-x_top(i-1))*((Pres_Foil_0_1(i)+Pres_Foil_0_1(i-1))/2)*tan(theta);
    F_Top_y(i-1) = -L*(x_top(i)-x_top(i-1))*((Pres_Foil_0_1(i)+Pres_Foil_0_1(i-1))/2);
    F_Top(i-1) = -sqrt((F_Top_x(i-1)^2+F_Top_y(i-1)^2));
    theta
end
F_Top_Total = sum(F_Top);

for i = 2:length(y_bot)
    theta = atan((y_bot(i)-y_bot(i-1))/(x_bot(i)-x_bot(i-1)));
    F_Bot_x(i-1) = L*(x_bot(i)-x_bot(i-1))*((Pres_Foil_0_1(i+10)+Pres_Foil_0_1(i-1+10))/2)*tan(theta);
    F_Bot_y(i-1) = -L*(x_bot(i)-x_bot(i-1))*((Pres_Foil_0_1(i+10)+Pres_Foil_0_1(i-1+10))/2);
    F_Bot(i-1) = sqrt((F_Bot_x(i-1)^2+F_Bot_y(i-1)^2));
    theta
end

F_Bot_y_tot = sum(F_Bot_y)
F_Top_y_tot = -sum(F_Top_y)



F_Bot_Total = sum(F_Bot)

F_Total = F_Bot_Total+F_Top_Total








