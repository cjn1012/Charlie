function acceleration = pointp(w1,w2)
%w1 must be in x and w1 must be in k

w1vector = [w1 0 0];
r  = [3*sqrt(3)/2 1.5 0];
w  = [w1 0 w2];
alpha = cross(w1vector,w);
vp = cross(w,r);


acceleration = cross(alpha,r) + cross(w,vp);
end











