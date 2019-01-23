% This script attempts to plot the velocity field of a 2D, fully devloped channel flow using the k-e model

% Constants
c_nu    = 0.09;
c_1      = 1.44;
c_2      = 1.92;
sigma_k = 1;
sigma_e = 1.3;
vis = 1.15*10^-5;

cells   = 100;
H       = 2;
rho     = 1;
u_tao   = 1;
delta_y = H/cells;

% Creating Array of u, k and epsilon
u   = ones(cells,1);
k   = ones(cells,1);
eps = ones(cells,1).*.01;

u(1)     = 0;
u(end)   = 0;
k(1)     = 0;
k(end)   = 0;
eps(1)   = 0;
eps(end) = 0;


% Iterating to get final values
tol = 1*10^-5;
err = 1;
for x = 1:cells*2
    for i = 2:cells-1
        u_old = u;
        k_old = k;
        eps_old = eps;
        eps_old(1) = eps_old(2);
        eps_old(end) = eps_old(end-1);
        v_t = c_nu*(k_old(i)^2/eps_old(i));
        u(i) = (delta_y/(2*(vis + v_t))) + (u_old(i+1)/2) + (u_old(i-1)/2);
        p_k = v_t * (((u_old(i+1)/2) + (u_old(i-1)/2))/delta_y);
        k(i) = (((p_k-eps_old(i))*delta_y)/(2*(vis+v_t/sigma_k)))+(k_old(i+1)/2) + (k_old(i-1)/2);
        eps(i) =    ((((eps_old(i)/k_old(i))*c_1*p_k)-c_2*(eps_old(i)^2)/k_old(i))*delta_y^2)/(2*(vis+v_t/sigma_e)) +(eps_old(i+1)/2) + (eps_old(i-1)/2);
    end
    err = abs(max(u)-max(u_old));
end

plot(u,linspace(0,2,cells))
xlabel('Flow Velocity [m/s]')
ylabel('Channel Height [m]')
