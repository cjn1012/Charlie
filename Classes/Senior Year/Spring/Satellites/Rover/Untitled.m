clear all
close all

L1 = linspace(1000,1400,900)
L2 = linspace(1000,1350,900)
H = 100
D = 1000

L1 = sqrt((L1.^2) + H.^2)
L2 = sqrt((L2.^2) + H.^2)

for i = 1:900
    A(i) = acosd(abs(((L2(i).^2)-(L1(i).^2)-(D.^2))/(-2.*L1(i).*D)))
end

if L1 > L2
    A = A+90
else
    A=A
end

for i = 1:900
    x(i) = L1(i)*cosd(A(i))-(D/2)
    y(i) = L1(i)*sind(A(i))
end

figure(1)
plot(x,y)
hold on
plot(-500,0,'o')
plot(500,0,'o')
title('Rover position over time')
xlabel('Longitude (m)')
ylabel('Latitude (m)')
xlim([-1000,1000])
ylim([0,2000])


