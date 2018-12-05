
V = 100
Air_D = .6
Q = 800000



NACA   = xlsread('xf-naca0015-il-500000.csv');
Alpha  = NACA(9:end,1);
Clv    = NACA(9:end,2);
Cdv    = NACA(9:end,3);
A = 10;
M_Equipment = 150;
M_Fuel_i = 150;
Fl_NoCl = .5*Air_D*V^2*A;
W_i = (1500+M_Equipment+M_Fuel_i)*9.8;

F_Drag = zeros(150,1);
Cd = zeros(150,1);
ranges = zeros(150,1);

for index = 1:150
    W = W_i - index*9.8;
    Cl = W/Fl_NoCl;
    for index2 = 1:150
        if Cl>Clv(index)
            Cd(index) = Cdv(index);
            x = Alpha(index);
        end
    end
    F_Drag(index) = .5 * Air_D * V^2 * A * Cd(index); % kgm3/s2 kgm2/s2
    
    ranges(index) = Q/F_Drag(index);
end

Range = sum(ranges);


