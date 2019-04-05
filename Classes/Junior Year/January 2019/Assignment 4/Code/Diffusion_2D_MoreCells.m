function [T] = Diffusion_2D_MoreCells(nx,ny)


%Specifying parameters
% Number of cells in each direction
%nx=100;               % Number of steps in space(x)
%ny=100;               % Number of steps in space(y) 

% Length of grid
L = 1.5;
H = 1.2;

% Time parameters for simulation
nt=100;              % Number of time steps 
dt=0.001;             % Width of each time step

% Defining the length of each cell sized linearly
dx=L/(nx-1);       % Width of space step(x)
dy=H/(ny-1);       % Width of space step(y)
x=0:dx:L;          % Range of x(0,1.5) and specifying the grid points
y=0:dy:H;          % Range of y(0,1.2) and specifying the grid points

% Assigning temperatures to the cells with T and Tn
T = zeros(nx,ny);       % Preallocating T
k = zeros(nx,ny);
% Assigning the diffustion coefficient/ viscocity numbers to use for the right sections of the grid of cells
vis_x_11_14 = 0.06;   % Diffusion coefficient k / viscocity
vis_y_05_07 = 0.06;   % Diffusion coefficient k / viscocity
vis_rem     = 22;     % Diffusion coefficient k / viscocity



%Initial Conditions of the cells

for i=1:nx
    for j=1:ny
        if (L/nx)*i > 1.1 && (L/nx)*i < 1.4                                   % Assigns the south boundary cells (y=0) to 16 degrees C
            k(i,j) = vis_x_11_14;
        elseif (H/ny)*j > 0.5 && (H/ny)*j < 0.7                                % Assigns the west boundary cells  (x=0) to 11 degrees C
            k(i,j) = vis_y_05_07;
        else                                        % Assigns the inside cells to 0 degrees C
            k(i,j) = vis_rem;
        end
    end
end

for i=1:nx
    for j=1:ny
        if j == 1                                  % Assigns the south boundary cells (y=0) to 16 degrees C
            T(i,j)=16;
        elseif i == nx                             % Assigns the west boundary cells  (x=0) to 11 degrees C
            T(i,j)= 11;
        elseif j == ny                             % Assigns the north boundary cells (y=1.2) to a varying degrees C
            T(i,j)=(20*(1+cos(pi*x(i)/L)));
        elseif i == 1                             % Assigns the east boundary cells  (x=1.5) to a varying degrees C
            T(i,j)=(12+20*(sin(pi*y(j)/H)));
        else                                        % Assigns the inside cells to 0 degrees C
            T(i,j)=0;  
        end
    end
end

% Solving for T with iterations
for x = 1:5000
    for i = 2:(nx-1)
        for j = 2:(ny-1)
            T(i,j) = (k(i+1,j)*T(i+1,j) + k(i-1,j)*T(i-1,j) + k(i,j+1)*T(i,j+1) + k(i,j-1)*T(i,j-1)) / (k(i+1,j) + k(i-1,j) + k(i,j+1) + k(i,j-1));
        end
    end
end



