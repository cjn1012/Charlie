%Charlie Nitschelm

% This script takes in parameters defining the configuration of the 2D diffusion problem
% and simulates its diffusion with various plots to help understand the situation
% through time as the temperature changes through time

% I USED ONE EQUATION FOR EVERY NODE! JUST MADE MATRIX VALUES FOR CONSTANTS WHERE THEY NEED TO BE
close all
clear all
%% Specifying parameters
f = 1;
Tavg(f) = 200; % Used to start while loop
f = 2;
Tavg(f) = 150;
L4 = 0; % Fin Length
while Tavg(f-1)-Tavg(f) > 10 % This constant determines how much change in T it will accept to keep going.
    f = f + 1; % Counter
    L4 = L4 + .01; % Increase fin length each iteration
    
    % Number of cells in each direction, 140 for nx and 25 for ny is done in around 2min.
    nx=145;              % Number of steps in space(x)
    ny=30;               % Number of steps in space(y)
    
    % Grid Dimensions
    L1 = 0.01; % Left Length
    L2 = 0.06; % Middle Length
    L3 = 0.01; % Right Length
    
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
    Index_1x =  ceil((L1/L)*nx);
    Index_2x =  floor(((L1+L2)/L)*nx);
    Index_3x =  floor(((L1+L2+L3)/L)*nx);
    Index_4x =  nx;
    
    Index_1y =  floor(((D1-F1)/D)*ny);
    Index_2y =  floor((D1/D)*ny);
    Index_3y =  ny;
    
    
    % Creating initial matricies for T and k in order to iterate over these to change T with each iteration to find steady state solution
    
    k = zeros(ny,nx);
    h = zeros(ny,nx);
    q = zeros(ny,nx);
    T = ones(ny,nx)*(30);
    q11 = 12000;
    
    for i=1:nx
        for j=1:ny
            if i > 1 && i<=Index_1x
                if j > 1 && j<Index_3y
                    k(j,i) = 260;
                end
            elseif i > Index_1x && i<=Index_2x
                if j > 1 && j<=Index_2y
                    k(j,i) = 260;
                end
            elseif i > Index_2x && i<=Index_3x
                if j > 1 && j<Index_3y
                    k(j,i) = 260;
                end
            elseif i > Index_3x && i<Index_4x
                if j >= Index_1y && j<=Index_2y
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
            elseif i >= Index_2x && i <= Index_3x+1
                if j == 1
                    h(j,i) = 180;
                elseif i == Index_3x+1
                    if j > 1 && j < Index_1y
                        h(j,i) = 180;
                    elseif j > Index_2y && j < Index_3y
                        h(j,i) = 180;
                    end
                end
            elseif i > Index_3x+1 && i<Index_4x
                if j == Index_1y-1
                    h(j,i) = 180;
                elseif j == Index_2y+1
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
            elseif i > Index_1x && i <= Index_2x
                if j == Index_2y
                    q(j,i) = q11;
                end
            elseif i == Index_2x+1
                if j >= Index_2y && j <= Index_3y
                    q(j,i) = q11;
                end
            end
            
        end
    end
    %
    T1 = T;
    T2 = T1+10;
    while abs(sum(sum(T2-T1))) > 1e-3
        T2 = T1;
        
        for j = 2:Index_1x
            for i = 2:(ny-1)
                T1(i,j) = GaussSolver(i,j,k,T1,h,q,dx,dy);
            end
        end
        
        for j = Index_1x+1:Index_2x
            for i = 2:Index_2y
                T1(i,j) = GaussSolver(i,j,k,T1,h,q,dx,dy);
            end
        end
        
        for j = Index_2x+1:Index_3x
            for i = 2:(ny-1)
                T1(i,j) = GaussSolver(i,j,k,T1,h,q,dx,dy);
            end
        end
        
        for j = Index_3x+1:Index_4x-1
            for i = Index_1y:Index_2y
                T1(i,j) = GaussSolver(i,j,k,T1,h,q,dx,dy);
            end
        end
        
    end
    Tavg(f) = mean(mean(T1));
end
Filepath = 'C:\Users\User\Desktop\Charlie\Classes\Junior Year\Spring 2019\Heat Transfer\Project 2\AsciiFile'; % make this your own!
AsciiMaker(Filepath,T1);
f1 = figure(1);
imagesc([0,L],[0,D],T1)
h = colorbar;
xlabel('Length (m)')
ylabel('Depth (m)')
ylabel(h, 'Temperature (°C)')
Filepath = 'C:\Users\User\Desktop\Charlie\Classes\Junior Year\Spring 2019\Heat Transfer\Project 2\TemperatureField'; % Make this your own!
saveas(f1, Filepath,'pdf');

