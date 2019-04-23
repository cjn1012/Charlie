%Problem 2

clear all
close all

J = 3;
k = 3;
bv = [0,3,6,12]; % Undamped, under,crital,over
tfinal = 15;
dt = .01;
timespan = linspace(0,tfinal,tfinal/dt);

b = bv(1);
sim('Problem2')
figure(1)
plot(tout,theta_0)
xlabel('Time (s)')
ylabel('Position ( )')
title('Undamped System')
grid on

b = bv(2);
sim('Problem2')
figure(2)
plot(tout,theta_0)
xlabel('Time (s)')
ylabel('Position ( )')
title('Underdamped System')
grid on

b = bv(3);
sim('Problem2')
figure(3)
plot(tout,theta_0)
xlabel('Time (s)')
ylabel('Position ( )')
title('Critically Damped System')
grid on

b = bv(4);
sim('Problem2')
figure(4)
plot(tout,theta_0)
xlabel('Time (s)')
ylabel('Position ( )')
title('Overdamped System')
grid on













