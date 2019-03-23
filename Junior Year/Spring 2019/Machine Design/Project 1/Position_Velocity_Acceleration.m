% Position Equations %

% Constants - Guessed

r_AB = 6;
r_BC = 4;
r_CD = 3;
r_DE = 5;


% Defined
r_EF = 1;
r_OA = 0.5;
r_FO = 3;
x_c  = 5;
y_c  = 2;

syms theta thetaa thetab phi
Loop_1x = r_OA*cos(theta) + r_AB*cos(thetaa) + r_BC*cos(thetab) - x_c == 0;
Loop_1y = r_OA*sin(theta) + r_AB*sin(thetaa) + r_BC*sin(thetab) - y_c == 0;
Loop_2x = r_CD*sin(thetab) + r_DE*cos(pi + phi) - r_EF +x_c == 0;
Loop_2y = r_CD*cos(thetab) + r_DE*sin(pi + phi) - r_FO + y_c == 0;

Solution = solve([Loop_1x,Loop_1y,Loop_2x,Loop_2y],[theta,thetaa,thetab,phi]);
ThetaSol = Solution.theta
ThetaaSol = Solution.thetaa
ThetabSol = Solution.thetab
PhiSol = Solution.phi
















