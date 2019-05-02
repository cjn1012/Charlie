
% This script takes in parameters defining the configuration of the 2D diffusion problem
% and simulates its diffusion with various plots to help understand the situation
% through time as the temperature changes through time

close all
clear all


%% Specifying parameters

% Number of cells in each direction
nx=50;               % Number of steps in space(x)
ny=50;               % Number of steps in space(y) 

% Grid Dimensions
L1 = 0.01; % Left Length
L2 = 0.06; % Middle Length
L3 = 0.01; % Right Length
L4 = 0.02; % Fin Length
D1 = 0.01; % Top Height
D2 = .015; % Bottom Height

L = L1+L2+L3+L4;
D = D1+D2;

% Fin Dimension
F1 = 0.005; % Fin Thickness


% Defining the length of each cell sized linearly
dx=L/(nx-1); % Width of space step(x)
dy=D/(ny-1);       % Width of space step(y)
x=0:dx:L;          % Range of x(0,1.5) and specifying the grid points
y=0:dy:D;          % Range of y(0,1.2) and specifying the grid points

% Assigning temperatures to the cells with T and Tn
Index_1x =  floor((L1/L)*nx);
Index_2x =  floor(((L1+L2)/L)*nx);
Index_3x =  floor(((L1+L2+L3)/L)*nx);
Index_4x =  nx;

Index_1y =  floor(((D1-F1)/D)*ny);
Index_2y =  floor((D1/D)*ny);
Index_3y =  ny;


%% Creating initial matricies for T and k in order to iterate over these to change T with each iteration to find steady state solution

k = zeros(ny,nx);
h = zeros(ny,nx);
q = zeros(ny,nx);
T = zeros(ny,nx)*(30+273);
T_00 = 30+273;
q11 = 12000;
for i=1:nx
    for j=1:ny
        if i >= 1 && i<=Index_1x
            if j >= 1 && j<=Index_3y                             
                T(j,i) = T_00;
            end
        elseif i > Index_1x && i<=Index_2x
            if j >= 1 && j<=Index_2y
                T(j,i) = T_00;
            end
        elseif i > Index_2x && i<=Index_3x
            if j >= 1 && j<=Index_3y  
                T(j,i) = T_00;
            end
        elseif i > Index_3x && i<=Index_4x
            if j >= Index_1y && j<=Index_2y  
                T(j,i) = T_00;
            end
        end
    end
end

for i=1:nx
    for j=1:ny
        if i > 1 && i<=Index_1x
            if j > 1 && j<Index_3y                             
                k(j,i) = 260;
            elseif j > Index_2y
                k(j,i) = 0;
            end
        elseif i > Index_1x && i<=Index_2x+1
            if j > 1 && j<=Index_2y  
                k(j,i) = 260;
            end
        elseif i > Index_2x+1 && i<=Index_3x
            if j > 1 && j<Index_3y-1  
                k(j,i) = 260;
            end
        elseif i > Index_3x && i<=Index_4x
            if j > Index_1y && j<Index_2y  
                k(j,i) = 260;
            end
        end
    end
end

for i=1:nx
    for j=1:ny
        if i > 1 && i<=Index_1x
            if j == 1                           
                h(j,i) = 180;
            end
        elseif i > Index_1x && i<Index_2x
            if j == 1 
                h(j,i) = 180;
            end
        elseif i >= Index_2x && i<=Index_3x
            if j == 1 
                h(j,i) = 180;
            elseif i == Index_3x
                if j > 1 && j < Index_1y
                    h(j,i) = 180;
                elseif j > Index_2y && j < Index_3y
                    h(j,i) = 180;
                end
            end
        elseif i > Index_3x && i<Index_4x
            if j == Index_1y
                h(j,i) = 180;
            elseif j == Index_2y
                h(j,i) = 180;
            end
        elseif i == Index_4x
            if j >= Index_1y && j <= Index_2y
                h(j,i) = 180;
            end
            
        end
    end
end

for i=1:nx
    for j=1:ny
        if i == Index_1x
            if j >= Index_2y && j <= Index_3y
                q(j,i) = q11;
            end
        elseif i > Index_1x && i < Index_2x
            if j == Index_2y
                q(j,i) = q11;
            end
        elseif i == Index_2x
            if j >= Index_2y && j <= Index_3y
                q(j,i) = q11;
            end
        end
        
    end
end

% a = dx;
% b = dy;
% c = k(i-1,j);
% d = k(i,j+1);
% e = k(i+1,j);
% f = k(i,j-1);
% h = T1(i-1,j);
% u = T1(i,j+1);
% z = T1(i+1,j);
% k = T1(i,j-1);
% q = h(i-1,j);
% r = h(i,j+1);
% s = h(i+1,j);
% t = h(i,j-1);
% m = q(i-1,j);
% n = q(i,j+1);
% o = q(i+1,j);
% p = q(i,j-1);
% 
% Solve = (a^2*b*h*q + a^2*b*z*s + a^2*b*m + a^2*b*o + a^2*c*h + a^2*e*z + a*b^2*u*r + a*b^2*k*t + a*b^2*n + a*b^2*p + b^2*d*u + ...
%     b^2*f*k)/(a^2*b*q + a^2*b*s + a^2*c + a^2*e + a*b^2*r + a*b^2*t + b^2*d + b^2*f)
%% Loop to iterate and calculate the T matrix 5000 times towards final convergence
T1 = T;

% J is going in the x direction, i is going down in the y direction


for x = 1:100
for j = 2:(nx-1)
    for i = 2:(ny-1)
        if j >= Index_1x && j <= Index_2x && i >= Index_2y
            flag = 1
            break
            
        elseif j >= Index_3x && i <= Index_1y && i >= Index_2y
            T1(i,j) = T1(i,j);
        end
        
        
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
        Solve(i,j) = (a^2*b*g*v + a^2*b*z*s + a^2*b*m + a^2*b*o + a^2*c*g + a^2*e*z + a*b^2*u*r + a*b^2*l*t + a*b^2*n + a*b^2*p + b^2*d*u + ...
            b^2*f*l)/(a^2*b*v + a^2*b*s + a^2*c + a^2*e + a*b^2*r + a*b^2*t + b^2*d + b^2*f);
        
        T1(i,j) = Solve(i,j);
        
        
    end
end

end






% for x = 1:10
%     
%     for j = 2:(nx-1)
%         for i = 2:(ny-1)
%             
%             Left Rectangle
%             a = dx;
%             b = dy;
%             c = k(i-1,j);
%             d = k(i,j+1);
%             e = k(i+1,j);
%             f = k(i,j-1);
%             g = T1(i-1,j);
%             u = T1(i,j+1);
%             z = T1(i+1,j);
%             l = T1(i,j-1);
%             v = h(i-1,j);
%             r = h(i,j+1);
%             s = h(i+1,j);
%             t = h(i,j-1);
%             m = q(i-1,j);
%             n = q(i,j+1);
%             o = q(i+1,j);
%             p = q(i,j-1);
%             Solve = (a^2*b*g*v + a^2*b*z*s + a^2*b*m + a^2*b*o + a^2*c*g + a^2*e*z + a*b^2*u*r + a*b^2*l*t + a*b^2*n + a*b^2*p + b^2*d*u + ...
%                     b^2*f*l)/(a^2*b*v + a^2*b*s + a^2*c + a^2*e + a*b^2*r + a*b^2*t + b^2*d + b^2*f);
%             if i == 1
%                 if j == 1
%                     T1(i,j) =Solve;
%                 elseif j > 1 && j < Index_3y
%                     T1(i,j) =Solve;
%                 elseif j == Index_3y
%                     T1(i,j) =Solve;
%                 end
%                 
%                 
%             elseif i > 1 && i < Index_1x
%                 if j == 1
%                     T1(i,j) =Solve;
%                 elseif j > 1 && j < Index_3y
%                     T1(i,j) = Solve;
%                 elseif j == Index_3y
%                     T1(i,j) = Solve;
%                 end
%                 
%                     
%             elseif i == Index_1x
%                 if j == 1
%                     T1(i,j) = Solve;
%                 elseif j > 1 && j < Index_2y
%                     T1(i,j) = Solve;
%                 elseif j == Index_2y
%                     T1(i,j) = Solve;
%                 elseif j > Index_2y && j < Index_3y
%                     T1(i,j) = Solve;
%                 elseif j == Index_3y
%                     T1(i,j) = Solve;
%                 end
%                 
% 
%                     
%             Middle Rectangle
%             
%             elseif i > Index_1x && i<Index_2x
%                 if j == 1
%                     T1(i,j) = Solve;
%                 elseif j > 1 && j < Index_2y
%                     T1(i,j) = Solve;
%                 elseif j == Index_2y
%                     T1(i,j) = Solve;
%                 end
%                 
%                     
%                     
%             Right Rectangle
%             
%             elseif i == Index_2x
%                 if j == 1
%                     T1(i,j) = Solve;
%                 elseif j > 1 && j < Index_2y
%                     T1(i,j) = Solve;
%                 elseif j >= Index_2y && j < Index_3y
%                     T1(i,j) = Solve;
%                 elseif j == Index_3y
%                     T1(i,j) = Solve;
%                 end
%                 
%                     
%             elseif i > Index_2x && i<Index_3x
%                 if j == 1
%                     T1(i,j) = Solve;
%                 elseif j > 1 && j < Index_3y
%                     T1(i,j) = Solve;
%                 elseif j == Index_3y
%                     T1(i,j) = Solve;
%                 end
%                 
%                     
%             elseif i == Index_3x
%                 if j == 1
%                     T1(i,j) = Solve;
%                 elseif j > 1 && j < Index_1y
%                     T1(i,j) = Solve;
%                 elseif j >= Index_1y && j < Index_2y
%                     T1(i,j) = Solve;        
%                 elseif j >= Index_2y && j < Index_3y
%                     T1(i,j) = Solve;
%                 elseif j == Index_3y
%                     T1(i,j) = Solve;
%                 end
%                 
%                     
%                 
%                 Fin
%             elseif i > Index_3x && i < Index_4x
%                 if j == Index_1y
%                     T1(i,j) = Solve;
%                 elseif j > Index_1y && j < Index_2y
%                     T1(i,j) = Solve;
%                 elseif j == Index_2y
%                     T1(i,j) = Solve;
%                 end
%                 
%                     
%             elseif i == Index_4x
%                 if j == Index_1y
%                     T1(i,j) = Solve;
%                 elseif j > Index_1y && j < Index_2y
%                     T1(i,j) = Solve;
%                 elseif j == Index_2y
%                     T1(i,j) = Solve;
%                 end
%                 
%                     
%             end
%         end
%     end
% end

figure(1)
imagesc(T1)
xlabel('Number of Cells (x-direction) [ ]')
ylabel('Number of Cells (y-direction) [ ]')
colorbar





