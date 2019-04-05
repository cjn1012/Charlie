% Position Equations %
degtorad = pi/180
% Constants - Guessed
theta = 5.71*degtorad;
thetaa = 7*pi/4;
thetab = pi/2;
phi = pi/8;



% Defined
r_EF = 2;
r_OA = 0.5;
r_FO = 3;
x_c  = 5;
y_c  = 2;

syms r_AB r_BC r_CD r_DE
Loop_1x = r_OA*cos(theta) + r_AB*cos(thetaa) + r_BC*cos(thetab) - x_c == 0;
Loop_1y = r_OA*sin(theta) + r_AB*sin(thetaa) + r_BC*sin(thetab) - y_c == 0;
Loop_2x = r_CD*sin(thetab) + r_DE*cos(pi + phi) - r_EF +x_c == 0;
Loop_2y = r_CD*cos(thetab) + r_DE*sin(pi + phi) - r_FO + y_c == 0;

Solution = solve([Loop_1x,Loop_1y,Loop_2x,Loop_2y],[r_AB,r_BC,r_CD,r_DE]);

r_AB = Solution.r_AB
r_BC = Solution.r_BC
r_CD = Solution.r_CD
r_DE = Solution.r_DE











