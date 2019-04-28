%%
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
L1 = 0.04; % Left Length
L2 = 0.06; % Middle Length
L3 = 0.04; % Right Length
L4 = 0.01; % Fin Length
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

Index_1y =  floor((D1/D)*ny);
Index_2y =  floor(((D1-F1)/D)*ny);
Index_3y =  ny;


%% Creating initial matricies for T and k in order to iterate over these to change T with each iteration to find steady state solution

%T = zeros(ny,nx);
k = ones(ny,nx)*10000000;

T_00 = 30+273;
%T_1  = 40+273;
for i=1:nx
    for j=1:ny
        if i >= 1 && i<=Index_1x
            if j >= 1 && j<=Index_3y                             
                T(j,i) = T_00;
            end
        elseif i > Index_1x && i<=Index_2x
            if j >= 1 && j<=Index_1y
                T(j,i) = T_00;
            end
        elseif i > Index_2x && i<=Index_3x
            if j >= 1 && j<=Index_3y  
                T(j,i) = T_00;
            end
        elseif i > Index_3x && i<=Index_4x
            if j >= Index_2y && j<=Index_1y  
                T(j,i) = T_00;
            end
        end
    end
end
          
for i=1:nx
    for j=1:ny
        if i > 1 && i<=Index_1x
            if j >= 1 && j<Index_3y                             
                k(j,i) = 260;
            end
        elseif i > Index_1x && i<=Index_2x
            if j >= 1 && j<=Index_1y  
                k(j,i) = 260;
            end
        elseif i > Index_2x && i<=Index_3x
            if j >= 1 && j<Index_3y  
                k(j,i) = 260;
            end
        elseif i > Index_3x && i<=Index_4x
            if j >= Index_2y && j<=Index_1y  
                k(j,i) = 260;
            end
        end
    end
end

%% Loop to iterate and calculate the T matrix 5000 times towards final convergence
T1 = T;
q = 12000; % W/m2
h = 180; % W/m2K

% J is going in the x direction, i is going down in the y direction

for x = 1:1000
    
    
for i = 2:(nx-1)
    for j = 2:(ny-1)


        % Left Rectangle
        if i > 1 && i < Index_1x
            if j == 1
                T1(j,i) = (T1(j,i) - (h*dx)/k(j,i));
                
            elseif j >= 2
                T1(j,i) = (k(j+1,i+1)*T1(j+1,i+1) + k(j+1,i-1)*T1(j+1,i-1) + k(j-1,i+1)*T1(j-1,i+1) + k(j-1,i-1)*T1(j-1,i-1)) / (k(j+1,i+1) + k(j+1,i-1) + k(j-1,i+1) + k(j-1,i-1));
            
            end
        elseif i == Index_1x
            if j >= Index_1y
                T1(j,i) = (T1(j,i) + (q*dy)/k(j,i));
            end
        
            
            
        %Middle Part
        elseif i >= Index_1x && i<Index_2x
            if j>=2 && j<Index_1y
                T1(j,i) = (k(j+1,i+1)*T1(j+1,i+1) + k(j+1,i-1)*T1(j+1,i-1) + k(j-1,i+1)*T1(j-1,i+1) + k(j-1,i-1)*T1(j-1,i-1)) / (k(j+1,i+1) + k(j+1,i-1) + k(j-1,i+1) + k(j-1,i-1));
            elseif j==1
                T1(j,i) = (T1(j,i) - (h*dx)/k(j,i));
            elseif j == Index_1y
                T1(j,i) = (T1(j,i) + (q*dx)/k(j,i));
            end
           
        % Right Rectangle  
        elseif i >= Index_2x && i<Index_3x
            if j>=2
                T1(j,i) = (k(j+1,i+1)*T1(j+1,i+1) + k(j+1,i-1)*T1(j+1,i-1) + k(j-1,i+1)*T1(j-1,i+1) + k(j-1,i-1)*T1(j-1,i-1)) / (k(j+1,i+1) + k(j+1,i-1) + k(j-1,i+1) + k(j-1,i-1));
            elseif j==1
                T1(j,i) = (T1(j,i) - (h*dx)/k(j,i));
            end
        elseif i == Index_2x
            if j >= Index_1y
                T1(j,i) = (T1(j,i) + (q*dy)/k(j,i));
            end
        elseif i == Index_3x
            if j >= 1 && j<=Index_2y
                T1(j,i) = (T1(j,i) - (h*dx)/k(j,i));
            elseif j >= Index_1y && j<=Index_3y
                T1(j,i) = (T1(j,i) - (h*dx)/k(j,i));
            end
            

        
        end
    end
end    
end
figure(1)
imagesc(T1)
xlabel('Number of Cells (x-direction) [ ]')
ylabel('Number of Cells (y-direction) [ ]')
colorbar





