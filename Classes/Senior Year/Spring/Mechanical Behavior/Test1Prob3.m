clear all
close all

G=120e9
b = .28e-9
v=0.3
L=1
x1 = linspace(-100,100,1000)
x2 = linspace(-10,10,1000)
j=0
for i = x1
    j=j+1
    f(j) = (G*b*L)/(2*pi*(1-v))*((3*i^2 + L^2)/((i^2 + L^2)^2))
end
k=0
for i = x2
    k=k+1
    h(k) = (G*b*L)/(2*pi*(1-v))*((3*i^2 + L^2)/((i^2 + L^2)^2))
end
figure(1)
plot(x1,f)
xlabel('Normalized Distance (x/L)')
ylabel('Force of Repulsion')

figure(2)
plot(x2,h)
xlabel('Normalized Distance (x/L)')
ylabel('Force of Repulsion')