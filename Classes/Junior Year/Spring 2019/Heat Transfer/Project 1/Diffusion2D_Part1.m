%%
% This script takes in parameters defining the configuration of the 2D diffusion problem
% and simulates its diffusion at steady state

close all
clear all


%% Specifying parameters

% Number of cells in each direction
nx=50;               % Number of steps in space(x)
ny=50;               % Number of steps in space(y) 

% Length of grid
L1 = 1.5;
L2 = 1.0;
D = 1.2;

% Assigning the diffustion coefficient/ viscocity numbers to use for
% section 1 and 2
k1 = 1;   % Diffusion coefficient k / viscocity
k2 = 0.5;   % Diffusion coefficient k / viscocity

% Defining the length of each cell sized linearly
dx=(L1+L2)/(nx-1);       % Width of space step(x)
dy=D/(ny-1);       % Width of space step(y)
x=0:dx:(L1+L2);          % Range of x(0,1.5) and specifying the grid points
y=0:dy:D;          % Range of y(0,1.2) and specifying the grid points

% Assigning temperatures to the cells with T and Tn
T = zeros(nx,ny);       % Preallocating T
k = zeros(nx,ny);



%% Creating initial matricies for T and k in order to iterate over these to change T with each iteration to find steady state solution


for i=1:nx
    for j=1:ny
        if j > (L2/L1)*nx  %((L1+L2)/nx)*j < L1         % Assigns the k values to their respective cell
            k(i,j) = k1;
        else                                      
            k(i,j) = k2;
        end
    end
end

T_0 = 20;
T_1 = 100;
for i=1:nx
    for j=1:ny
        if j == 1                                  % Assigns the east boundary cells (y=0) to 16 degrees C
            T(i,j) = T_0;
        elseif i == nx                             % Assigns the south boundary cells  (x=0) to 11 degrees C
            T(i,j) = T_0;
        elseif i == 1
            if j > (L2/L1)*nx                % Assigns the north boundary cells (y=1.2) to a varying degrees C
                T(i,j) = T_0;
            else
                T(i,j) = T_1;
            end
        elseif j == ny                              % Assigns the west boundary cells  (x=1.5) to a varying degrees C
            T(i,j)= T_0;
        else                                       % Assigns the inside cells to 0 degrees C
            T(i,j)=0;  
        end
    end
end

%% Loop to iterate and calculate the T matrix 5000 times towards final convergence
T2 = T;

for x = 1:1000
    for i = 2:(nx-1)
        for j = 2:(ny-1)
            T2(i,j) = (k(i+1,j)*T2(i+1,j) + k(i-1,j)*T2(i-1,j) + k(i,j+1)*T2(i,j+1) + k(i,j-1)*T2(i,j-1)) / (k(i+1,j) + k(i-1,j) + k(i,j+1) + k(i,j-1));
        end
    end
end




% 
% for int = 1:5000
%     for i = 2:(nx-1)
%         for j = 2:(ny-1)
%             T2(i,j) = (k(i+1,j)*T(i+1,j) + k(i-1,j)*T(i-1,j) + k(i,j+1)*T(i,j+1) + k(i,j-1)*T(i,j-1)) / (k(i+1,j) + k(i-1,j) + k(i,j+1) + k(i,j-1)); % Calculate the new T matrix
%         end
%     end
% end    
%     
%     % These if statements with for loops helps calculate flux between iteration 199 and 200 by calculating difference of T and then apply k to it 
%     if int == 199
%         for i = 2:(nx-1)
%             for j = 2:(ny-1)
%                 TotalFlux_x(i,j) = (k(i+1,j)*T(i+1,j) + k(i-1,j)*T(i-1,j)) / (k(i+1,j) + k(i-1,j));
%                 TotalFlux_y(i,j) = (k(i,j+1)*T(i,j+1) + k(i,j-1)*T(i,j-1)) / (k(i,j+1) + k(i,j-1));
%                 TotalFlux(i,j) = (k(i+1,j)*T(i+1,j) + k(i-1,j)*T(i-1,j) + k(i,j+1)*T(i,j+1) + k(i,j-1)*T(i,j-1)) / (k(i+1,j) + k(i-1,j) + k(i,j+1) + k(i,j-1));
%             end
%         end
%     end
%     if int == 200
%         for i = 2:(nx-1)
%             for j = 2:(ny-1)
%                 TotalFlux_x(i,j) = k(i,j)*(TotalFlux_x(i,j)-((k(i+1,j)*T(i+1,j) + k(i-1,j)*T(i-1,j)) / (k(i+1,j) + k(i-1,j))));
%                 TotalFlux_y(i,j) = k(i,j)*(TotalFlux_y(i,j)-((k(i,j+1)*T(i,j+1) + k(i,j-1)*T(i,j-1)) / (k(i,j+1) + k(i,j-1))));
%                 TotalFlux(i,j) = k(i,j)*(TotalFlux(i,j)-((k(i+1,j)*T(i+1,j) + k(i-1,j)*T(i-1,j) + k(i,j+1)*T(i,j+1) + k(i,j-1)*T(i,j-1)) / (k(i+1,j) + k(i-1,j) + k(i,j+1) + k(i,j-1))));
%             end
%         end
%     end
%     T = T2;
% end

% %% Plotting Data T
% 
% 
% % Mesh Convergence Plots
% figure(1)
% imagesc(T)
% xlabel('Number of Cells (x-direction) [ ]')
% ylabel('Number of Cells (y-direction) [ ]')
% colorbar
% 
% figure(2)
% imagesc(Diffusion_2D_MoreCells(100,100))
% xlabel('Number of Cells (x-direction) [ ]')
% ylabel('Number of Cells (y-direction) [ ]')
% colorbar
% 
% figure(3)
% imagesc(Diffusion_2D_MoreCells(500,500))
% xlabel('Number of Cells (x-direction) [ ]')
% ylabel('Number of Cells (y-direction) [ ]')
% colorbar
% 
% % Contour Plot of Temperature
% figure(4)
% contour(T,'ShowText','on')
% xlabel('Number of Cells (x-direction) [ ]')
% ylabel('Number of Cells (y-direction) [ ]')
% 
% 
% % 3D Color Plot of Temperature
% figure(5)
% contour3(T)
% xlabel('Number of Cells (x-direction) [ ]')
% ylabel('Number of Cells (y-direction) [ ]')
% zlabel('Temperature [°C]')
% 
% 
% % Vector Plot of the Heat Flux
% figure(6)
% [x,y] = meshgrid(x(1:49),y(1:49));
% quiver(x,y,TotalFlux_x,TotalFlux_y);
% xlabel('Distance [m]')
% ylabel('Distance [m]')
% 
% % Heat Flux 3D Color Plot
% figure(7)
% surf(x,y,TotalFlux)
% xlabel('Distance [m]')
% ylabel('Distance [m]')
% zlabel('Total Flux [W/m^2]')
% % 2 Additional Plots
% figure(8)
% imagesc(TotalFlux)
% xlabel('Number of Cells (x-direction) [ ]')
% ylabel('Number of Cells (y-direction) [ ]')
% colorbar
% 
% 
% figure(9)
% contour3(Diffusion_2D_NoNorthHT(100,100))
% xlabel('Number of Cells (x-direction) [ ]')
% ylabel('Number of Cells (y-direction) [ ]')
% zlabel('Temperature [°C]')
% 
% figure(13)
% imagesc(Diffusion_2D_NoNorthHT(100,100))
% xlabel('Number of Cells (x-direction) [ ]')
% ylabel('Number of Cells (y-direction) [ ]')
% colorbar
% 
% 
% % For Modified BC, 3D Contour Plot of Temperature
% figure(10)
% imagesc(Diffusion_2D_Iterations(50))
% xlabel('Number of Cells (x-direction) [ ]')
% ylabel('Number of Cells (y-direction) [ ]')
% colorbar
% 
% figure(11)
% imagesc(Diffusion_2D_Iterations(500))
% xlabel('Number of Cells (x-direction) [ ]')
% ylabel('Number of Cells (y-direction) [ ]')
% colorbar
% 
% figure(12)
% imagesc(Diffusion_2D_Iterations(5000))
% xlabel('Number of Cells (x-direction) [ ]')
% ylabel('Number of Cells (y-direction) [ ]')
% colorbar
% 
% 
% 
% 
