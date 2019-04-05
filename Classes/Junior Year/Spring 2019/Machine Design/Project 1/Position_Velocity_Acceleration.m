clear all
close all

%%
% Position Equations %


% Constants - Guessed

r_AB = 5.160;
r_CD = 1.12;
r_BC = 1.12;
r_DE = 2.95;

% Defined
r_EF = 1;
r_OA = 0.5;
r_FO = 3;
x_c  = 5;
y_c  = 2;

% Position Equations for each Point

i=1;
for theta = 0:360
    xA(i) = r_OA*cosd(theta);
    yA(i) = r_OA*sind(theta);
    i = i+1;
end

i=1;
for theta_1 = 203:1:333
    xB(i) = 5 + r_BC*cosd(theta_1);
    yB(i) = 2 + r_BC*sind(theta_1);
    i = i+1;
end

i=1;
for theta_2 = 203:1:333
    xD(i) = 5 - r_CD*cosd(theta_2);
    yD(i) = 2 - r_CD*sind(theta_2);
    i = i+1;
end

i=1;
for d = 0:1
    xE(i) = 1.5 + d;
    yE(i) = 3;
    i = i+1;
end

figure(1)
plot(xA,yA,xB,yB,xD,yD,xE,yE)
xlim([-1,7])
ylim([-1,5.5])




%% Graph 2

i = 1;
for theta = 0:.01:360
    xA(i) = cosd(theta);
    xE1(i) = .0074+cosd(theta-185.7);
    i = i+1;
end

figure
hold on
theta = linspace(0,360,36001);
plot(theta,xA,theta,xE1)
ylim([-2,10])
xlim([0,370])
title('X Position vs Crank Angle')
legend('X Position of Point A', 'X Position of point E')
hold off 


%% Graph 3

i = 1;
for theta = 0:.01:360
    xA(i) = -sind(theta);
    xE2(i) = .0074-sind(theta-185.7);
    i = i+1;
end

figure
hold on
theta = linspace(0,360,36001);
plot(theta,xA,theta,xE2)
ylim([-2,10])
xlim([0,370])
title('X Velocity vs Crank Angle')
legend('X Velocity of Point A', 'X Velocity of point E')
hold off 

%% Graph 4

i = 1;
for theta = 0:.01:360
    xA3(i) = -cosd(theta);
    xE3(i) = .0074-cosd(theta-185.7);
    aOA(i) = sqrt((0.25*cosd(theta))^2 + (0.25*sind(theta))^2);
    
    i = i+1;
end


figure
hold on
theta = linspace(0,360,36001);
plot(theta,xA3,theta,xE3)
ylim([-2,10])
xlim([0,370])
title('X Acceleration vs Crank Angle')
legend('X Acceleration of Point A', 'X Acceleration of point E')
hold off 

%%






















% syms theta thetaa thetab phi
% Loop_1x = r_OA*cos(theta) + r_AB*cos(thetaa) + r_BC*cos(thetab) - x_c == 0;
% Loop_1y = r_OA*sin(theta) + r_AB*sin(thetaa) + r_BC*sin(thetab) - y_c == 0;
% Loop_2x = r_CD*sin(thetab) + r_DE*cos(pi + phi) - r_EF +x_c == 0;
% Loop_2y = r_CD*cos(thetab) + r_DE*sin(pi + phi) - r_FO + y_c == 0;
% 
% Solution = solve([Loop_1x,Loop_1y,Loop_2x,Loop_2y],[theta,thetaa,thetab,phi]);
% ThetaSol = Solution.theta
% ThetaaSol = Solution.thetaa
% ThetabSol = Solution.thetab
% PhiSol = Solution.phi
















