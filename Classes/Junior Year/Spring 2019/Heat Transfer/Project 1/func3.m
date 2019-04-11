function [A , B, cc, rr] = func3(L1,L2, D, T_0, T_1,k1,k2);

% Number of cells in each direction
res = 100
nx=res*(L1+L2);               % Number of steps in space(x)
ny=res*D;               % Number of steps in space(y) 


% Assigning temperatures to the cells with T and Tn
T = zeros(nx,ny);       % Preallocating T
k = zeros(nx,ny);


for i=1:nx
    for j=1:ny
        if j > (L2/L1)*nx    % Assigns the k values to their respective cell
            k(i,j) = k1;
        else                                      
            k(i,j) = k2;
        end
    end
end

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

for x = 1:5000
    for i = 2:(nx-1)
        for j = 2:(ny-1)
            T2(i,j) = (k(i+1,j)*T2(i+1,j) + k(i-1,j)*T2(i-1,j) + k(i,j+1)*T2(i,j+1) + k(i,j-1)*T2(i,j-1)) / (k(i+1,j) + k(i-1,j) + k(i,j+1) + k(i,j-1));
        end
    end
end








