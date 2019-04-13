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

AoA_Foil_0_T  = [AoA_Foil_0(1),AoA_Foil_0(3:12),AoA_Foil_0(12)];
AoA_Foil_9_T  = [AoA_Foil_0(1),AoA_Foil_0(3:12),AoA_Foil_0(12)];
AoA_Foil_18_T = [AoA_Foil_0(1),AoA_Foil_0(3:12),AoA_Foil_0(12)];

AoA_Foil_0_B  = [AoA_Foil_0(1),AoA_Foil_0(13:22),AoA_Foil_0(22)];
AoA_Foil_9_B  = [AoA_Foil_0(1),AoA_Foil_0(13:22),AoA_Foil_0(22)];
AoA_Foil_18_B = [AoA_Foil_0(1),AoA_Foil_0(13:22),AoA_Foil_0(22)];


AoA_Foil_0_T_1  = AoA_Foil_0_T-(AoA_Foil_0(end));
AoA_Foil_9_T_1  = AoA_Foil_9_T-(AoA_Foil_9(end));
AoA_Foil_18_T_1 = AoA_Foil_18_T-(AoA_Foil_18(end));

AoA_Foil_0_B_1  = AoA_Foil_0_B-(AoA_Foil_0(end));
AoA_Foil_9_B_1  = AoA_Foil_9_B-(AoA_Foil_9(end));
AoA_Foil_18_B_1 = AoA_Foil_18_B-(AoA_Foil_18(end));


Pres_Foil_0_1_T =  -Den_Water*Gravity.*AoA_Foil_0_T_1;
Pres_Foil_9_1_T =  -Den_Water*Gravity.*AoA_Foil_9_T_1;
Pres_Foil_18_1_T = -Den_Water*Gravity.*AoA_Foil_18_T_1;
Pres_Foil_0_1_B =  -Den_Water*Gravity.*AoA_Foil_0_B_1;
Pres_Foil_9_1_B =  -Den_Water*Gravity.*AoA_Foil_9_B_1;
Pres_Foil_18_1_B = -Den_Water*Gravity.*AoA_Foil_18_B_1;

Pres_Foil_0_T =  Pres_Foil_0_1_T(1:12);
Pres_Foil_9_T =  Pres_Foil_9_1_T(1:12);
Pres_Foil_18_T = Pres_Foil_18_1_T(1:12);
Pres_Foil_0_B =  Pres_Foil_0_1_B(1:12);
Pres_Foil_9_B =  Pres_Foil_9_1_B(1:12);
Pres_Foil_18_B = Pres_Foil_18_1_B(1:12);


c = .15; % Chorde Length
t = .12;% Thickness to chord
L = .30; % Span

x_top = [0,0.76, 3.81, 11.43, 19.05, 38, 62, 80.77, 101.35, 121.92, 137.16,150]./1000;
x_bot = [0,1.52,7.62,15.24,22.86,41.15,59.44,77.73,96.02,114.3,129.54,150]./1000; 

y_top =5*c*t.*( (.2969.*sqrt(x_top/c))-(.1260.*x_top/c)-(.3516.*(x_top/c).^2)+(.2843.*(x_top/c).^3)-(.1015.*(x_top/c).^4));
y_bot =-5*c*t.*( (.2969.*sqrt(x_bot/c))-(.1260.*x_bot/c)-(.3516.*(x_bot/c).^2)+(.2843.*(x_bot/c).^3)-(.1015.*(x_bot/c).^4));

plot(x_top,y_top,x_bot,y_bot)

% Pres_Foil_0_T = ones(1, length(y_top));
% Pres_Foil_0_B = ones(1, length(y_top));


% Top Stuff
for i = 1:length(y_top)-1
    theta = atan((y_top(i+1)-y_top(i))/(x_top(i+1)-x_top(i)));
    Drag_Top_0(i)  = L*(x_top(i+1)-x_top(i))*((Pres_Foil_0_T(i)+Pres_Foil_0_T(i+1))/2)*tan(theta);
    Drag_Top_9(i)  = L*(x_top(i+1)-x_top(i))*((Pres_Foil_9_T(i)+Pres_Foil_9_T(i+1))/2)*tan(theta);
    Drag_Top_18(i)  = L*(x_top(i+1)-x_top(i))*((Pres_Foil_18_T(i)+Pres_Foil_18_T(i+1))/2)*tan(theta);
    Lift_Top_0(i)  = -L*(x_top(i+1)-x_top(i))*((Pres_Foil_0_T(i)+Pres_Foil_0_T(i+1))/2);
    Lift_Top_9(i)  = -L*(x_top(i+1)-x_top(i))*((Pres_Foil_9_T(i)+Pres_Foil_9_T(i+1))/2);
    Lift_Top_18(i)  = -L*(x_top(i+1)-x_top(i))*((Pres_Foil_18_T(i)+Pres_Foil_18_T(i+1))/2);

    b_theta = atan((y_bot(i+1)-y_bot(i))/(x_bot(i+1)-x_bot(i)));
    Drag_Bot_0(i)  = -L*(x_bot(i+1)-x_bot(i))*((Pres_Foil_0_B(i)+Pres_Foil_0_B(i+1))/2)*tan(b_theta);
    Drag_Bot_9(i)  = -L*(x_bot(i+1)-x_bot(i))*((Pres_Foil_9_B(i)+Pres_Foil_9_B(i+1))/2)*tan(b_theta);
    Drag_Bot_18(i) = -L*(x_bot(i+1)-x_bot(i))*((Pres_Foil_18_B(i)+Pres_Foil_18_B(i+1))/2)*tan(b_theta);
    Lift_Bot_0(i)  = L*(x_bot(i+1)-x_bot(i))*((Pres_Foil_0_B(i)+Pres_Foil_0_B(i+1))/2);
    Lift_Bot_9(i)  = L*(x_bot(i+1)-x_bot(i))*((Pres_Foil_9_B(i)+Pres_Foil_9_B(i+1))/2);
    Lift_Bot_18(i) = L*(x_bot(i+1)-x_bot(i))*((Pres_Foil_18_B(i)+Pres_Foil_18_B(i+1))/2);
end

d_top_tot_0  = sum(Drag_Top_0)
d_bot_tot_0  = sum(Drag_Bot_0)

l_top_tot_0  = sum(Lift_Top_0)
l_bot_tot_0  = sum(Lift_Bot_0)










