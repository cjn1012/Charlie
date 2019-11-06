clear all
close all


% Part b


% P-Control

num1 = [168.4]
den1 = [1,189.6]
figure(1)
sys = tf(num1,den1)
rlocus(sys)

% I-Control
num2 = [34048]
den2 = [1,189.6,0]
figure(2)
sys2 = tf(num2,den2)
rlocus(sys2)

% PI-Control
num3 = [168.8,666.7*168.7]
den3 = [1,189.6,0]
figure(3)
sys3 = tf(num3,den3)
rlocus(sys3)