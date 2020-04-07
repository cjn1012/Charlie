clear all
close all


ix =100;
iy = 100;
iz = 1000;

sim('model')

plot(tout,u(:,5))
