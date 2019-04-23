%% Problem 3

clear all

% Constants
J = 3;
K = 3;
B = 3;

% Time stuff
t_i = 0;
t_f = 20;
dt = 0.01;
step = t_f/dt;


% State space vars
A = [0,1;-K/J,-B/J];
B = [0;K/J];
C = [1,0;0,1];
D = [0;0];
x_0 = [-1;0];
u = 1;
x_dot = A*x_0 + B*u;

for t = 1:step
    x_0 = x_dot*dt + x_0;
    x_dot = A*x_0 + B*u;
    y(t) = x_0(1);
end

time = linspace(0,step*dt,step);

figure(5)
plot(time,y)
title('Underdamped System: J=3,K=3,B=3')
grid on
xlabel('Time (s)')
ylabel('Position')

