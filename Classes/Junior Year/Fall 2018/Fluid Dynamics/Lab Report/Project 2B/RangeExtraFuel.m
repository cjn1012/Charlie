function [Range, efficiency3] = RangeExtraFuel(V,Air_D,Q,M)
Air_Data = xlsread('air_data1.xls');
Densities = Air_Data(1:43,10);
Pressures = Air_Data(1:43,8);
Temperatures = Air_Data(1:43,7);
Efficiency2 = .5019*.9;
Q = Q*Efficiency2;
NACA   = xlsread('xf-naca0015-il-1000000.csv');
Alpha  = NACA(12:end,1);
Clv    = NACA(12:end,2);
Cdv    = NACA(12:end,3);
A = 25;
M_Equipment = 150;
M_Fuel_i  = M;
Fl_NoCl = .5*Air_D*V^2*A;
W_i = (1500+M_Equipment+M_Fuel_i)*9.78;
ranges = zeros(M,1);
efficiency3 = Efficiency(20,1500,interp1(Densities,Pressures,Air_D),interp1(Densities,Temperatures,Air_D)+273);
for index = 1:1:M
    W = W_i - index*9.78;
    Cl = W/Fl_NoCl;
    Cd = interp1(Clv,Cdv,Cl);
    F_Drag = .5 * Air_D * V^2 * A * Cd; % kgm3/s2 kgm2/s2
    ranges(index) = Q/F_Drag;
    
end
Range = nansum(ranges)/1000;
end
