function [x] = func2(A, B)


B = B';   %transpose(B)
x = A\B;

end
