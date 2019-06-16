L = 1
k = 240
T1 = 600
T2 = 300
r1 = .1
r2 = 1
r = linspace(r1,r2,100);
for i = 1:length(r)
    T(i) = ((T1-T2)/(log(r1/r2)))*log(r(i)/r2)+T2;
end

plot(r, T)

