function solve = Solver(i,j,k,T1,h,q,dx,dy)
a = dx;
b = dy;
c = k(i-1,j);
d = k(i,j+1);
e = k(i+1,j);
f = k(i,j-1);
g = T1(i-1,j);
u = T1(i,j+1);
z = T1(i+1,j);
l = T1(i,j-1);
v = h(i-1,j);
r = h(i,j+1);
s = h(i+1,j);
t = h(i,j-1);
m = q(i-1,j);
n = q(i,j+1);
o = q(i+1,j);
p = q(i,j-1);
solve = (a^2*b*g*v + a^2*b*z*s + a^2*b*m + a^2*b*o + a^2*c*g + a^2*e*z + a*b^2*u*r + a*b^2*l*t + a*b^2*n + a*b^2*p + b^2*d*u + ...
    b^2*f*l)/(a^2*b*v + a^2*b*s + a^2*c + a^2*e + a*b^2*r + a*b^2*t + b^2*d + b^2*f);


end
