close all

K = 52
J = .01
B = .13


plot(theta)
xlim([0,1])
ylabel('Theta')
hold on
yyaxis right
plot(theta_dot)
ylabel('Theta Dot')

hold off
figure(2)
Ri = 10000
Rf = 1000000
c = .000010

E = tf([-1],[c*Ri,Ri/Rf,0])
bode(E)
