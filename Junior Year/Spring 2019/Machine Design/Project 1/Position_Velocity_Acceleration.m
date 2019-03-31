% Position Equations %
clear all
close all
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
xlim([0,7])
ylim([0,7])

%%%%%%%%%%%%%%%%%%%%%







i = 1;
for theta = 0:.01:360
    xA(i) = sind(theta);
    xE1(i) = sind(theta-185.7);
    i = i+1;
end
theta = linspace(0,360,36001);
plot(theta,xA,theta,xE1)

xlim([0,360])





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
















