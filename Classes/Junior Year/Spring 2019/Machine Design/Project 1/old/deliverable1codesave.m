clear all; close all;

%% For Loops for Calling Functions
theta = (0:1:360).*(pi/180);
for i= 1:length(theta)
    x0 = [5.7*pi/180,120*pi/180,0*pi/180,2.5];        
    options=optimset('Display','iter');  
    x(i,:) = fsolve(@(x) myfun(x,theta(i)),x0,options); 
end
for i= 1:length(theta)
    x0 = [5.7*pi/180,120*pi/180,0*pi/180,2.5];        
    options=optimset('Display','iter');  
    y(i,:) = fsolve(@(y) myfun1(y,theta(i),x(i,:)),x0,options);         
end

for i= 1:length(theta)
    x0 = [5.7*pi/180,120*pi/180,0*pi/180,2.5];        
    options=optimset('Display','iter');  
    z(i,:) = fsolve(@(z) myfun2(z,theta(i),x(i,:),y(i,:)),x0,options);         
end

r_AB = 5.160;
r_CD = 1.12;
r_BC = 1.12;
r_DE = 2.95;

% Defined
r_EF = 1;
r_OA = 0.5;
r_FO = 3;
x_c  = 5;
y_c  = 2;

theta_dot = 50;


%Loop #1
R_OAx=(r_OA/2)*cos(theta);
R_OAy=(r_OA/2)*sin(theta);
R_ABx=(r_AB/2)*cos(x(:,1));
R_ABy=(r_AB/2)*sin(x(:,1));
R_BCx=(r_BC/2)*cos(x(:,2));
R_BCy=(r_BC/2)*sin(x(:,2));
R_COx=(x_c/2);
R_COy=(y_c/2);

%Loop #2
R_CDx=(r_CD/2)*sin(x(:,2));
R_CDy=(r_CD/2)*cos(x(:,2));
R_DEx=(r_DE/2)*cos(pi+x(:,3));
R_DEy=(r_DE/2)*sin(pi+x(:,3));
R_EFx=(1/2);
R_EFy=(0/2);
R_FOx=(0/2);
R_FOy=(1/2);


%%

%position

function F = myfun(x,theta)
r_AB = 5.160;
r_CD = 1.12;
r_BC = 1.12;
r_DE = 2.95;

% Defined
r_EF = 1;
r_OA = 0.5;
r_FO = 3;
x_c  = 5;
y_c  = 2;

theta_dot = 50;

    F = [r_OA*cos(theta)+r_AB*cos(x(:,1))+r_AB*cos(x(:,2))-x_c;
    r_OA*sin(theta)+r_AB*sin(x(:,1))+r_AB*sin(x(:,2))-y_c;
    x_c+r_CD*sin(x(:,2))+r_DE*cos(pi+x(:,3))-x(:,4);
    y_c+r_CD*cos(x(:,2))+r_DE*sin(pi+x(:,3))-r_FO];
   
    
    end 

%velocity

    function F = myfun1(y,theta,x)
r_AB = 5.160;
r_CD = 1.12;
r_BC = 1.12;
r_DE = 2.95;

% Defined
r_EF = 1;
r_OA = 0.5;
r_FO = 3;
x_c  = 5;
y_c  = 2;

theta_dot = 50;
F = [-r_OA*theta_dot*sin(theta)-r_AB*y(1)*sin(x(:,1))-r_BC*y(2)*sin(x(:,2));
    r_OA*theta_dot*cos(theta)+r_AB*y(1)*cos(x(:,1))+r_AB*y(2)*cos(x(:,2));
    r_CD*y(2)*cos(x(:,2))-r_DE*y(3)*sin(pi+x(:,3));
    -r_CD*y(2)*sin(x(:,2))+r_DE*y(3)*cos(pi+x(:,3))-y(4)];
    
    end 

%%% acceleration

function F = myfun2(z,theta,x,y)
r_AB = 5.160;
r_CD = 1.12;
r_BC = 1.12;
r_DE = 2.95;

% Defined
r_EF = 1;
r_OA = 0.5;
r_FO = 3;
x_c  = 5;
y_c  = 2;

theta_dot = 50;
F = [-r_OA*(theta_dot^2)*cos(theta)-r_AB*z(1)*sin(y(:,1))-r_AB*(y(:,1)^2)*cos(y(:,1))-r_BC*z(2)*sin(y(:,2))-r_BC*(y(:,2)^2)*cos(y(:,2));
    -r_OA*(theta_dot^2)*sin(theta)+r_AB*z(1)*cos(y(:,1))-r_AB*(y(:,1)^2)*sin(y(:,1))+r_BC*z(2)*cos(y(:,2))-r_BC*(y(:,2)^2)*cos(y(:,2));
    r_CD*z(2)*cos(y(:,2))-r_CD*(y(:,2)^2)*sin(y(:,2))-r_DE*z(3)*sin(pi+x(:,3))-r_DE*(y(:,3)^2)*cos(pi+x(:,3));
    -r_CD*z(2)*sin(y(:,2))-r_CD*(y(:,2)^2)*cos(y(:,2))+r_DE*z(3)*cos(pi+x(:,3))-r_DE*(y(:,3)^2)*sin(pi+x(:,3))-z(4)];

end



% function F = myfun3(w,x_2,x_4,z_4, R_O1_x, R_O1_y,R_1B_x,R_1B_y,...
%         R_B2_x, R_B2_y,R_A2_x, R_A2_y,R_23_x, R_23_y, R_34_x,...
%         R_34_y,R_32_x, R_32_y)
% 
%     
% %k = 1;  %lbf/in
% %x_spring = x_4;
% %x_spring = x
% %F_spring = k*(-x_4-1.3125);  %%% Spring Force %%%
%    
% 
% 
% F = [w(1) + w(3);
%     w(2) + w(4);
%          
%     w(5) + w(2)*R_O1_x - w(1)*R_O1_y + w(4)*R_1B_x - w(3)*R_1B_y;
% 
%     -w(3) + w(6)*cos(pi/2 - x_2);
%     -w(4) - w(6)*sin(pi/2 - x_2);
% 
%     w(7) - w(6)*cos(pi/2 - x_2) + w(9);
%     w(8) + w(6)*sin(pi/2 - x_2) + w(10);
% 
%     w(6)*cos(pi/2 - x_2)*R_B2_y - w(6)*sin(pi/2 - x_2)*R_B2_x...
%     + w(10)*R_A2_x - w(9)*R_A2_y + w(8)*R_23_x - w(7)*R_23_y;
%     
% 
%     w(11) - w(7);
%     w(12) - w(8);
% 
%     w(12)*R_34_x - w(11)*R_34_y - w(8)*R_32_x + w(7)*R_32_y;
% 
%     1*(-x_4-1.3125) - w(11) - (27/9650)*z_4;
%     w(13) - w(12);]
%  
% 
%  
% end
