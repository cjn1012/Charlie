%Problem 2

clear all
close all
M = 2
J = 2
k = 2
b = 0

tfinal = 15
dt = .1

timespan = linspace(1,tfinal,tfinal/dt)
sim('Problem2')

figure(1)
plot(timespan,theta_0)





















