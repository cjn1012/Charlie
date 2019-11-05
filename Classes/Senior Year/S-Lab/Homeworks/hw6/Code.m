clear all
close all

%Constants
Ktach = .003
R = 4.2
J = .0004
B = .002
Ka = -1
Kt = 6.602
Ke = .047
Ri = 100000
Rf = 150000
C1 = .000000033
C2 = .00000001

% Part b

% p-control

num1 = [Rf*Kt*Ktach/(B*R)]
den1 = [Ri*J/B,Ri]
sys = tf(num1,den1)



rlocus(sys)
