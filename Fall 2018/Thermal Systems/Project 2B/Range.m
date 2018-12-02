function Range = Range(V)


Efficiency = .6016;
Q_Comb = 930; %J/kg

NACA = xlsread('xf-naca0015-il-500000.csv');
Alpha = NACA(12:end,1);
Clv    = NACA(12:end,2);
Cdv    = NACA(12:end,3);

Air_D = .6;
A = 1;
M_Equipment =  150;
M_Fuel = 100;
W = 1500+M_Equipment+M_Fuel;


Fl_NoCl = .5*Air_D*V^2*A;

Cl = W/Fl_NoCl;

for index = 1:150
    if Cl>Clv(index)
        Cd = Cdv(index);
        x = Alpha(index)
        break
    end
end

F_Drag = .5 * Air_D * V^2 * A * Cd; % kgm3/s2 kgm2/s2
W_d = Q_Comb * Efficiency; %J/kg
%kgm

Range =  (W_d / F_Drag)*1000*(M_Fuel*(((W-M_Fuel)/W)+1))










end







