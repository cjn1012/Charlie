%THIS FUNCTION FINDS THE LARGEST/SMALLEST CELL RATIO GIVEN:
%INPUT length of curve l, cell expansion ratio cer, smallest cell size delta
%OUTPUT largest/smallest cell ratio Rratio

function R=Rratio(l,r,delta)

R=l*(r-1)/delta/r+1/r;
