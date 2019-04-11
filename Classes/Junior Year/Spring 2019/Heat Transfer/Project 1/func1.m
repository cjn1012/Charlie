
function [A , B, cc, rr] = func1(L1,L2, D, T_0, T_1,k1,k2,ResX,ResY)
L = L1+L2;

cc = ResX;
rr = ResY; 
dx = (L1+L2)/(cc-1);
dy = (D)/(rr-1);
    
ccL1 = round((L1/L)*ResX);

NodeType = zeros(ResX,ResY);

NodeType(1,1) = 1; % Top-Left Node
NodeType(1, 2:ccL1) = 2; % Right Top Nodes
NodeType(1, ccL1+1:cc-1) = 3; % Left Top Nodes
NodeType(1, cc) = 4; % Top-Right Node

NodeType(2:rr-1, 1) = 5; % Middle-Left Node
NodeType(2:rr-1, 2:ccL1) = 6; % Right Middle Nodes
NodeType(2:rr-1, ccL1+1:cc-1) = 7; % Left Middle Nodes
NodeType(2:rr-1, cc) = 8; % Middle-Right Node

NodeType(rr, 1) = 9; % Bottom-Left Node
NodeType(rr, 2:ccL1) = 10; % Right Bottom Nodes
NodeType(rr, ccL1+1:cc-1) = 11; % Left Bottom Nodes
NodeType(rr, cc) = 12; % Bottom-Right Node




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
                
            case 3  %% Node 10 , Top Middle right
                A(n,n) = -( ((2*dx)/dy) + (dy/dx));
                A(n,n+1) = (dy/(2*dx)); %% Right
                A(n,n-1) = (dy/(2*dx)); %% Right
                A(n,n+cc) = (dx/dy);    %% Down node 2
                B(n) = -(dx/dy)*T_0;  %% Known value(s)
                
            case 4  %% Node 3 , Top Right
                A(n,n) = -( (dx/dy) + (dy/dx) );
                A(n,n-1) = (dy/(2*dx));     %% Left
                A(n,n+cc) = (dx/(2*dy));    %% Down node 
                B(n) = -(dx/dy)*T_0; %% Known value(s)   
                
            case 5  %% Node 4 , Left Middle
                A(n,n) = -( (dx/dy) + ((2*dy)/dx) );
                A(n,n+1) = (dy/dx);         %% Right
                A(n,n+cc) = (dx/(2*dy));    %% Down node 
                A(n,n-cc) = (dx/(2*dy));    %% Upper node
                B(n) = -(dy/dx)*T_0;          %% Known value(s)  
                
            case 6  %% Node 5 , Middle
                A(n,n) = -( ((2*dx)/dy) + ((2*dy)/dx) );
                A(n,n+1) = (dy/dx);     %% Right
                A(n,n-1) = (dy/dx);     %% Left
                A(n,n+cc) = (dx/dy);    %% Down node 
                A(n,n-cc) = (dx/dy);    %% Upper node
                B(n) = 0;               %% Known value(s)  
                
            case 7  %% Node 5 , Middle right side
                A(n,n) = -( ((2*dx)/dy) + ((2*dy)/dx) );
                A(n,n+1) = (dy/dx);     %% Right
                A(n,n-1) = (dy/dx);     %% Left
                A(n,n+cc) = (dx/dy);    %% Down node 
                A(n,n-cc) = (dx/dy);    %% Upper node
                B(n) = 0;               %% Known value(s) 
                
            case 8  %% Node 5 , Middle-Right Node
                A(n,n) = -( ((2*dy)/dx) + (dx/dy) );
                A(n,n-1) = (dy/dx);         %% Left
                A(n,n+cc) = (dx/(2*dy));    %% Down node 
                A(n,n-cc) = (dx/(2*dy));    %% Upper node
                B(n) = -(dy/dx)*T_0;               %% Known value(s) 
                
            case 9 
                A(n,n) = -( (dy/dx) + (dx/dy) );
                A(n,n+1) = (dy/(2*dx));     %% Right
                A(n,n-cc) = (dx/(2*dy));    %% Upper node
                B(n) = -( (dy/(2*dx)) + (dx/(2*dy)) )*T_0;    %% Known value(s) 
                
            case 10 
                A(n,n) = -( ((2*dx)/dy) + (dy/dx) );
                A(n,n+1) = (dy/(2*dx));     %% Right
                A(n,n-1) = (dy/(2*dx));     %% Left
                A(n,n-cc) = (dx/dy);    %% Upper node
                B(n) = -(dx/dy)*T_0;    %% Known value(s) 
                
            case 11 
                A(n,n) = -( ((2*dx)/dy) + (dy/dx) );
                A(n,n+1) = (dy/(2*dx));     %% Right
                A(n,n-1) = (dy/(2*dx));     %% Left
                A(n,n-cc) = (dx/dy);   
                B(n) = -(dx/dy)*T_0;   
                
            case 12
                A(n,n) = -( (dx/dy) + (dy/dx) );
                A(n,n-1) = (dy/(2*dx));      
                A(n,n-cc) = (dx/(2*dy));      
                B(n) = -( (dy/(2*dx)) + (dx/(2*dy)) )*T_0; 
                
        end
    end
    
end


end