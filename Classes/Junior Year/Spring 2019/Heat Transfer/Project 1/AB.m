function [x] = AB(A, B)

B = B';   %transpose(B)
x = A\B;

end
