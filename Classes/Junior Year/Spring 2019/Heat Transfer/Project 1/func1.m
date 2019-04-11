
function [A , B, cc, rr] = func1(L1,L2, D, T_0, T_1,k1,k2,ResX,ResY);

cc = ResX*(L1+L2); 
rr = ResY*(D); 

dx = 1/(rr-1);
dy = 1/(cc-1);

NodeType = zeros(rr,cc);

NodeType(1,1) = 1;
NodeType(1, 2:L1) = 2;
NodeType(1, L1+1:cc-1) = 10;
NodeType(1, cc) = 3;

NodeType(2:rr-1, 1) = 4;
NodeType(2:rr-1, 2:L1) = 5;
NodeType(2:rr-1, L1+1:cc-1) = 11;
NodeType(2:rr-1, cc) = 6;

NodeType(rr, 1) = 7;
NodeType(rr, 2:L1) = 8;
NodeType(rr, L1+1:cc-1) = 12;
NodeType(rr, cc) = 9;




for i = 1:rr
    for j = 1:cc
       n = j + ((i-1)*cc);
        
       switch NodeType(i,j)
            
            case 1  %% Node 1 , Top Left
                A(n,n) = -((dy/dx) + (dx/dy));  %% Node #
                A(n,n+1) = (dy/(2*dx));     %% Right
                A(n,n+cc) = (dx/(2*dy));    %% Down node 2
                B(n) = -( (dy/(2*dx))*T_0 + (dx/(2*dy))*T_1);  %% Known values
                
            case 2  %% Node 2 , Top left middle stuff
                A(n,n) = -( ((2*dx)/dy) + (dy/dx));
                A(n,n+1) = (dy/(2*dx)); %% Right
                A(n,n-1) = (dy/(2*dx)); %% Right
                A(n,n+cc) = (dx/dy);    %% Down node 2
                B(n) = -(dx/dy)*T_1;  %% Known value(s)
                
            case 3  %% Node 3 , Top Right
                A(n,n) = -( (dx/dy) + (dy/dx) );
                A(n,n-1) = (dy/(2*dx));     %% Left
                A(n,n+cc) = (dx/(2*dy));    %% Down node 
                B(n) = -( (dy/(2*dx))*T_0 + (dx/(2*dy))*T_1); %% Known value(s)   
                
            case 4  %% Node 4 , Left Middle
                A(n,n) = -( (dx/dy) + ((2*dy)/dx) );
                A(n,n+1) = (dy/dx);         %% Right
                A(n,n+cc) = (dx/(2*dy));    %% Down node 
                A(n,n-cc) = (dx/(2*dy));    %% Upper node
                B(n) = -(dy/dx)*T_0;          %% Known value(s)  
                
            case 5  %% Node 5 , Middle
                A(n,n) = -( ((2*dx)/dy) + ((2*dy)/dx) );
                A(n,n+1) = (dy/dx);     %% Right
                A(n,n-1) = (dy/dx);     %% Left
                A(n,n+cc) = (dx/dy);    %% Down node 
                A(n,n-cc) = (dx/dy);    %% Upper node
                B(n) = 0;               %% Known value(s)  
                
            case 6  %% Node 5 , Middle Right
                A(n,n) = -( ((2*dy)/dx) + (dx/dy) );
                A(n,n-1) = (dy/dx);         %% Left
                A(n,n+cc) = (dx/(2*dy));    %% Down node 
                A(n,n-cc) = (dx/(2*dy));    %% Upper node
                B(n) = -(dy/dx)*T_0;               %% Known value(s) 
                
            case 7 
                A(n,n) = -( (dy/dx) + (dx/dy) );
                A(n,n+1) = (dy/(2*dx));     %% Right
                A(n,n-cc) = (dx/(2*dy));    %% Upper node
                B(n) = -( (dy/(2*dx)) + (dx/(2*dy)) )*T_0;    %% Known value(s) 
                
            case 8 
                A(n,n) = -( ((2*dx)/dy) + (dy/dx) );
                A(n,n+1) = (dy/(2*dx));     %% Right
                A(n,n-1) = (dy/(2*dx));     %% Left
                A(n,n-cc) = (dx/dy);    %% Upper node
                B(n) = -(dx/dy)*T_0;    %% Known value(s) 
                
            case 9
                A(n,n) = -( (dx/dy) + (dy/dx) );
                A(n,n-1) = (dy/(2*dx));        %% Left
                A(n,n-cc) = (dx/(2*dy));       %% Upper node
                B(n) = -( (dy/(2*dx)) + (dx/(2*dy)) )*T_0;    %% Known value(s) 
                
            case 10  %% Node 10 , Top Middle right
                A(n,n) = -( ((2*dx)/dy) + (dy/dx));
                A(n,n+1) = (dy/(2*dx)); %% Right
                A(n,n-1) = (dy/(2*dx)); %% Right
                A(n,n+cc) = (dx/dy);    %% Down node 2
                B(n) = -(dx/dy)*T_0;  %% Known value(s)
                            
            case 11  %% Node 5 , Middle right side
                A(n,n) = -( ((2*dx)/dy) + ((2*dy)/dx) );
                A(n,n+1) = (dy/dx);     %% Right
                A(n,n-1) = (dy/dx);     %% Left
                A(n,n+cc) = (dx/dy);    %% Down node 
                A(n,n-cc) = (dx/dy);    %% Upper node
                B(n) = 0;               %% Known value(s) 
                
            case 12 
                A(n,n) = -( ((2*dx)/dy) + (dy/dx) );
                A(n,n+1) = (dy/(2*dx));     %% Right
                A(n,n-1) = (dy/(2*dx));     %% Left
                A(n,n-cc) = (dx/dy);    %% Upper node
                B(n) = -(dx/dy)*T_0;    %% Known value(s) 
        end
    end
    
end

cc = L1+L2
rr = D
end