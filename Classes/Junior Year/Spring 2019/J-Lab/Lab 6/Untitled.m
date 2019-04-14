%% 
close all; clc; clear all;

Ti          = 21.2+273.15;                   % [C]
Tf          = 22.6+273.15; 
Pi          = 102268;                        % [Pa]
Pf          = 102235;    
R           = 8.314;                         % [J/molK]
MMass       = 28.97;                         % [g/m3]
Di_air      = (MMass*Pi)/(R*Ti*1000);        % [kg/m3]
Df_air      = (MMass*Pf)/(R*Tf*1000); 
D_w         = 998.2;
g           = 9.81;
chord       = 150/1000;
L        = 300/1000;
Planf_A     = chord*L;
t           = .12; 

% Part 2

AoA_0 = [1.59, 3.2, 2.95, 3.54, 3.85, 3.57, 3.79, 3.65, 3.52,...+
    3.05, 3.31, 3.16, 2.98, 3.55, 3.65, 3.44, 3.21, 3.5, 3.4,...+
    3.31, 2.91, 3.14, .49]./39.37;              % [m]
    
AoA_9  = [1.6, 3.2, 7, 5.9, 5.25, 4.4, 4.2, 3.8, 3.55,...+
    3.05, 3.25, 3.15, 1.6, 2.2, 2.6, 2.65, 2.7, 3.05, 3.1,...+
    3.05, 2.7, 2.9, .45]./39.37;                % [m]

AoA_18 = [1.6, 3.25, 4.05, 3.95, 4.00, 3.7, 4.1, 4.15, 4.2,...+
    3.65, 4.0, 3.85, 1.5, 1.9, 2.25, 2.35, 2.5, 2.9, 3.0, 3.1,...+
    2.85, 3.25, 0.5]./39.37;                    % [m]

AirfoilPress = Pi - g*D_w.*AoA_0;

% Part 3
% Change in pressure used to determine fluid velocity 
dP_0    = D_w*g*(AoA_0(2) - AoA_0(1));
dP_9    = D_w*g*(AoA_9(2) - AoA_9(1));
dP_18   = D_w*g*(AoA_18(2) - AoA_18(1));

Vel_0   = sqrt((2/Di_air)*dP_0);
Vel_9   = sqrt((2/Di_air)*dP_9);
Vel_18   = sqrt((2/Di_air)*dP_18);
Vel = [Vel_0, Vel_9, Vel_18];
Vel = mean(Vel);

% Part 4

AoA = -9:3:21;
D_force = [1.5, 1.1, 0.8, 0.75, 0.85, 1.1, 1.5, 1.95, 2.5, 5.1, 5.7];
L_force = [-9.3, -6.3, -3, 0.4, 3.3, 6.3, 9, 11.5, 13.5, 11.1, 10.8];

CL = L_force./(0.5*Di_air*Vel^2*Planf_A);
CD = D_force./(0.5*Di_air*Vel^2*Planf_A);

p = polyfit(AoA, L_force, 1);
L_pred = p(1)*AoA + p(2);
zeroL = -p(2)/p(1);                 % Don't Need to adjust


% Part 5
figure(1)
plot(AoA, CD, AoA, CL)
xlabel('Angle of Attack (AoA) [\circ]')
ylabel('Coefficient of Drag/Lift [ ]')
legend('C_{D}', 'C_{L}', 'location', 'southeast')

% Part 6
Coeff_Ratio = CL./CD;

figure(2)
plot(AoA, Coeff_Ratio)
xlabel('Angle of Attack (AoA) [\circ]')
ylabel('C_{D}/C_{L} [ ]')

%% Part 7
close all; clc;

D_w         = 998.2;
g           = 9.81;
chord       = 150/1000;
L        = 300/1000;
Planf_A     = chord*L;
t           = .12; 

xt = [0, 0.76, 3.81, 11.43, 19.05, 38, 62, 80.77, 101.35, 121.92, 137.16, 150]./1000;
xb = [0, 1.52, 7.62, 15.24, 22.86, 41.15, 59.44, 77.73, 96.02, 114.3, 129.54, 150]./1000;

a = xt./chord;
b = xb./chord;

yt = 5.*t.*chord.*((0.2969.*sqrt(a)) + (-.1260.*(a))...+
    + (-.3516.*(a).^2) + (.2843.*(a).^3) + (-.1015.*(a).^4));
yb = -5.*t.*chord.*((0.2969.*sqrt(b)) + (-.1260.*(b))...+
    + (-.3516.*(b).^2) + (.2843.*(b).^3) + (-.1015.*(b).^4));

figure(3)
plot(xt, yt, xb, yb)
xlabel('[meters]'); ylabel('[meters]'); title('Airfoil')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0 DEGREE ANGLE OF ATTACK %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta_t    = zeros(1,length(yt) - 1);
theta_b    = zeros(1,length(yb) - 1);
DragTop0 = zeros(1,length(yt) - 1);
LiftTop0 = zeros(1,length(yt) - 1);
DragBot0 = zeros(1,length(yb) - 1);
LiftBot0 = zeros(1,length(yb) - 1);

Force_xt_9 = zeros(1,length(yt) - 1);
Force_yt_9 = zeros(1,length(yt) - 1);
Force_xb_9 = zeros(1,length(yb) - 1);
Force_yb_9 = zeros(1,length(yb) - 1);

Force_xt_18 = zeros(1,length(yt) - 1);
Force_yt_18 = zeros(1,length(yt) - 1);
Force_xb_18 = zeros(1,length(yb) - 1);
Force_yb_18 = zeros(1,length(yb) - 1);


height_0t = AoA_0(end) - [AoA_0(1), AoA_0(3:12), AoA_0(12)];
height_0b = AoA_0(end) - [AoA_0(1), AoA_0(13:22), AoA_0(22)]; 

Pres_0t = height_0t.*g*D_w;
Pres_0b = height_0b.*g*D_w;
% Pres_0t = ones(1, length(height_0t));
% Pres_0b = ones(1, length(height_0b));

height_9t = AoA_9(end) - [AoA_9(1), AoA_9(3:12), AoA_9(12)];
height_9b = AoA_9(end) - [AoA_9(1), AoA_9(13:22), AoA_9(22)]; 

Pres_9t = height_9t.*g*D_w;
Pres_9b = height_9b.*g*D_w;
% Pres_9t = ones(1, length(height_9t));
% Pres_9b = ones(1, length(height_9t));

height_18t = AoA_18(end) - [AoA_18(1), AoA_18(3:12), AoA_18(12)];
height_18b = AoA_18(end) - [AoA_18(1), AoA_18(13:22), AoA_18(22)]; 

Pres_18t = height_18t.*g*D_w;
Pres_18b = height_18b.*g*D_w;

for i = 1:length(yt)-1
    % Top Calculations
    theta_t(i)     =  atand( (yt(i+1) - yt(i)) / (xt(i+1) - xt(i)) );
    DragTop0(i)    =  L*(xt(i+1) - xt(i) ) * ( (Pres_0t(i+1)  + Pres_0t(i) ) /2)*tand(theta_t(i));
    Force_xt_9(i)  =  L*(xt(i+1) - xt(i) ) * ( (Pres_9t(i+1)  + Pres_9t(i) ) /2)*tand(theta_t(i));
    Force_xt_18(i) =  L*(xt(i+1) - xt(i) ) * ( (Pres_18t(i+1) + Pres_18t(i)) /2)*tand(theta_t(i));
    LiftTop0(i)    = -L*(xt(i+1) - xt(i) ) * ( (Pres_0t(i+1)  + Pres_0t(i) ) /2);
    Force_yt_9(i)  = -L*(xt(i+1) - xt(i) ) * ( (Pres_9t(i+1)  + Pres_9t(i) ) /2);
    Force_yt_18(i) = -L*(xt(i+1) - xt(i) ) * ( (Pres_18t(i+1) + Pres_18t(i)) /2);
    
    
    % Bottom Calculations
    theta_b(i)     =  atand((yb(i+1) - yb(i))/(xb(i+1) - xb(i)));
    DragBot0(i)    = - L*(xb(i+1) - xb(i) ) * ( (Pres_0b(i+1)  + Pres_0b(i))/2)*tand(theta_b(i));
    Force_xb_9(i)  = - L*(xb(i+1) - xb(i) ) * ( (Pres_9b(i+1)  + Pres_9b(i))/2)*tand(theta_b(i));
    Force_xb_18(i) = - L*(xb(i+1) - xb(i) ) * ( (Pres_18b(i+1) + Pres_18b(i))/2)*tand(theta_b(i));
    LiftBot0(i)    =  L*(xb(i+1) - xb(i) ) * ( (Pres_0b(i+1)  + Pres_0b(i))/2);
    Force_yb_9(i)  =  L*(xb(i+1) - xb(i) ) * ( (Pres_9b(i+1)  + Pres_9b(i))/2);    
    Force_yb_18(i) =  L*(xb(i+1) - xb(i) ) * ( (Pres_18b(i+1) + Pres_18b(i))/2);
end

totalDrag_0 = sum(DragTop0 + DragBot0)
totalLift_0 = sum(LiftTop0 + LiftBot0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 9 DEGREE ANGLE OF ATTACK %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = 9;

Force_Dt9 =  Force_xt_9.*cosd(alpha) + Force_yt_9.*sind(alpha);
Force_Lt9 = -Force_xt_9.*sind(alpha) + Force_yt_9.*cosd(alpha);

Force_Db9 =  Force_xb_9.*cosd(alpha) + Force_yb_9.*sind(alpha);
Force_Lb9 =  -Force_xb_9.*sind(alpha) + Force_yb_9.*cosd(alpha);

totalDrag_9 = sum(Force_Dt9 + Force_Db9) 
totalLift_9 = sum(Force_Lt9 + Force_Lb9)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 18 DEGREE ANGLE OF ATTACK %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = 18;

Force_Dt18 =  Force_xt_18.*cosd(alpha) + Force_yt_18.*sind(alpha);
Force_Lt18 = -Force_xt_18.*sind(alpha) + Force_yt_18.*cosd(alpha);

Force_Db18 = Force_xb_18.*cosd(alpha) + Force_yb_18.*sind(alpha);
Force_Lb18 = -Force_xb_18.*sind(alpha) + Force_yb_18.*cosd(alpha);

totalDrag_18 = sum(Force_Dt18 + Force_Db18)
totalLift_18 = sum(Force_Lt18 + Force_Lb18)

Fd_0 = [0; totalDrag_0];
Fd_9 = [9; totalDrag_9];
Fd_18 = [18; totalDrag_18];

T = table(Fd_0, Fd_9, Fd_18);

%% Part 8

nu = 1.516e-5;

Re = (Vel*chord)/nu;
skinFric = 0.0742/(Re^(1/5));
F_skin = skinFric*(0.5*Di_air*Vel^2*Planf_A);

%% Part 9 
ar = L/chord; 
eff_fact = 0.7;

CDi = CL.^2/(pi*ar*eff_fact);
Ind_drag = CDi.*(0.5*Di_air*Vel^2*Planf_A); 

ind_0 = [0; Ind_drag(4)];
ind_9 = [9; Ind_drag(7)];
ind_18 = [18; Ind_drag(10)];

T2 = table(ind_0, ind_9, ind_18);
T3 = table(AoA', Ind_drag');

figure(4)
plot(AoA, Ind_drag, '-o')
xlabel('Angle of Attack')
ylabel('Induced Drag, F_{D,i} [N]')

%% Part 11


T4 = table([0 9 18]',[Ind_drag(4), Ind_drag(7), Ind_drag(10)]',[totalDrag_0 totalDrag_9 totalDrag_18]',[D_force(4), D_force(7), D_force(10)]');

TOTAL_0 = F_skin + Ind_drag(4) + totalDrag_0;
TOTAL_9 = F_skin + Ind_drag(7) + totalDrag_9;
TOTAL_18 = F_skin + Ind_drag(10) + totalDrag_18;
TOTAL = [TOTAL_0 TOTAL_9 TOTAL_18];

figure(5)
plot([AoA(1) AoA(end)], [F_skin F_skin], AoA, Ind_drag,...+
    [0 9 18], [totalDrag_0 totalDrag_9 totalDrag_18],...+
    AoA, D_force, [0 9 18], TOTAL);
xlabel('Angle of Attack [\circ]')
ylabel('Drag Force [N]')
legend('Est. Skin Drag', 'Est. Induced Drag', 'Est. Form Drag',...+
    'Recorded Drag Force', 'Est. Total Drag', 'location', 'northwest');

%% Part 12
u_L = 0.01; 
u_A = 0.25*10^-6;  
u_Man = 0.05;
u_D = 0.05;
 
d_L = 1/(u_A*D_w*g*(u_Man(1) - u_Man(2)));
d_A = -Lift(9)/((u_A^2)*D_w*g*(u_Man(1) - u_Man(2)));
d_Man = -Lift(9)/(u_A*D_w*g*(u_Man(1) - u_Man(2))^2);
d_D = -Lift(9)/(u_A*(D_w)^2*g*(u_Man(1) - u_Man(2)));
 
uncertainty = sqrt((d_L*u_L)^2+(d_Man*u_Man)^2+(d_D*u_D)^2+(d_A*u_A)^2);












