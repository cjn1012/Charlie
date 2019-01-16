%THIS CODE PRINTS OUT THE COORDINATES FOR A SPLINE (EXCLUDING THE INITIAL
%AND FINAL VERTEX).
%INPUT: in the input section you have to provide:
%the x, y and z coordinates of the points, or alternatively: 
%the x, y, or z coordinates of the start and end vertices, 
%the number of desired  and
%modify the equations to provide the relationships between x, y, and z
%coordinates
%OUTPUT: coordinates of points placed in paranthesis

%% INPUT OPTION 1
%please enter x, y and z coordinates of the points below
%make sure that the sizes of x, y, and z are the same
N=1:1000;
x=1;    %e.g. x=N/(N+1);            %the reason why you would divide by N+1 is to avoid having point of x=1 at the end if this is you vertex coordinate
y=1;    %e.g. y=x.^2                %the dot means that you multiply x element-wise rather than multiiplying it as a matrix
z=1;    %e.g. z=2*ones(length(x))   %this makes z the same size as x, and all the values are 2    

%% INPUT OPTION 2
if length(x)==1     %this just checks whether option 1 was used
    %HOW MANY POINTS DO YOU WANT TO GENERATE?
    N=1000;  
    %PROVIDE x, y or z COORDINATE OF FIRST VERTEX
    x1=0;                  
    %PROVIDE COORDINATE OF SECOND VERTEX
    x2=2;                  
    
    x=x1+(x2-x1)*(1:N)/(N+1);   %this will distribute the points evenly in x direction
    
    %PROVIDE RELATIONSHIP BETWEEN x, y, and z
    y=x.^2+x*4;                 %the dot means that the operation is performed element-wise
    z=1.5*ones(length(x));      %fixes z to 1.5 for all points
end
    
%% THIS WRITES A TEXT FILE WITH THE COORDINATES
f=fopen('spline','w');

for j=1:length(x)
    fprintf(f, '(%f %f %f) \n', x(j), y(j), z(j));
end

fclose(f);