V = 200;
Eff_Comb = 1;
Air_D = .6;
Q = 720;


Efficiency = .5129;
Q_Comb = Eff_Comb*Q; %J/kg

NACA   = xlsread('xf-naca0015-il-500000.csv');
Alpha  = NACA(12:end,1);
Clv    = NACA(12:end,2);
Cdv    = NACA(12:end,3);
A = 10;
M_Equipment =  150;
M_Fuel_i = 150;
Fl_NoCl = .5*Air_D*V^2*A;
W_i = 1500+M_Equipment+M_Fuel_i;



ranges = zeros(150,1);

for index = 1:150
    W = W_i - index;
    Cl = W/Fl_NoCl;
    if Cl>Clv(index)
        Cd = Cdv(index);
        x = Alpha(index);
    end
    F_Drag = .5 * Air_D * V^2 * A * Cd; % kgm3/s2 kgm2/s2
    ranges(index) = 42000000/F_Drag;
end

Range = sum(ranges);














