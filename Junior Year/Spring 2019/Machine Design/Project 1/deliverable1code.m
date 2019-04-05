clear all; close all;

%% For Loops for Calling Functions
theta = (0:1:360).*(pi/180);
for i= 1:length(theta)
    x0 = [5.7*pi/180,120*pi/180,175*pi/180,2.5];        
    options=optimset('Display','iter');  
    pos(i,:) = fsolve(@(pos) myfun(pos,theta(i)),x0,options); 
end
for i= 1:length(theta)
    x0 = [5.7*pi/180,120*pi/180,175*pi/180,2.5];        
    options=optimset('Display','iter');  
    vel(i,:) = fsolve(@(vel) myfun1(vel,theta(i),pos(i,:)),x0,options);         
end

for i= 1:length(theta)
    x0 = [5.7*pi/180,120*pi/180,175*pi/180,2.5];        
    options=optimset('Display','iter');  
    acc(i,:) = fsolve(@(acc) myfun2(acc,theta(i),pos(i,:),vel(i,:)),x0,options);         
end

r_AB = 5.160;
r_CD = 1.12;
r_BC = 1.12;
r_DE = 2.95;

% Defined
r_EF = 2.5;
r_OA = 0.5;
r_FO = 3;
x_c  = 5;
y_c  = 2;

theta_dot = 50;


%Loop #1
R_OAx=(r_OA/2)*cos(theta);
R_OAy=(r_OA/2)*sin(theta);
R_ABx=(r_AB/2)*cos(pos(:,1));
R_ABy=(r_AB/2)*sin(pos(:,1));
R_BCx=(r_BC/2)*cos(pos(:,2));
R_BCy=(r_BC/2)*sin(pos(:,2));
R_COx=(x_c/2);
R_COy=(y_c/2);

%Loop #2
R_CDx=(r_CD/2)*sin(pos(:,2));
R_CDy=(r_CD/2)*cos(pos(:,2));
R_DEx=(r_DE/2)*cos(pi+pos(:,3));
R_DEy=(r_DE/2)*sin(pi+pos(:,3));
R_EFx=(1/2);
R_EFy=(0/2);
R_FOx=(0/2);
R_FOy=(1/2);

%% Dynamics Graphs
%% Graphs

r_EFy = r_FO; 
yEF(1:361) = r_EFy;
figure;
plot(pos(:,4),yEF, 'Linewidth', 1)
hold on
plot(.5*cos(theta),.5*sin(theta),'Linewidth', 1)
plot(1.12*cos(pos(:,2))+x_c , 1.12*sin(pos(:,2))+y_c, 'Linewidth', 1)
plot(1.12*cos(pos(:,2)+pi)+x_c , 1.12*sin(pos(:,2)+pi)+y_c, 'Linewidth', 1)
xlabel(' X [inches] ')
ylabel(' Y [inches] ')
legend('r EF','r OA','point D','point B','Location','northwest')
axis([-1 6 -1 3.5])

figure;
title('X pos of (A&E) vs Crank Angle') 
hold on 
plot(theta*(180/pi),.5*cos(theta),theta*(180/pi),pos(:,4))
legend('A','E')
hold off

figure 
hold on
title('X vel of (A&E) vs Crank Angle')
plot(theta*(180/pi),-.5*sin(theta)*theta_dot)
plot(theta*(180/pi),vel(:,4))
legend('A','E')
hold off

figure 
hold on
title('X acc of (A&E) vs Crank Angle')
plot(theta*(180/pi),-.5*cos(theta))
plot(theta*(180/pi),acc(:,4))
legend('A','E')
hold off

% figure;
% %plot(phi1,a_cm_OB)
% plot(phi1,a_cm_AC)
% hold on
% plot(phi1,a_cm_CD)
% plot(phi1, a_cm_block)
% xlabel('Crank Angle [rad]')
% ylabel('Magnitude of Acceleration [in/s^2]')
% legend('AC','CD','Block')


%%

%position

function F = myfun(pos,theta)
r_AB = 5.160;
r_CD = 1.12;
r_BC = 1.12;
r_DE = 2.95;

% Defined
r_EF = 2.5;
r_OA = 0.5;
r_FO = 3;
x_c  = 5;
y_c  = 2;

theta_dot = 50;

    F = [r_OA*cos(theta)+r_AB*cos(pos(:,1))+r_BC*cos(pos(:,2))-x_c;
    r_OA*sin(theta)+r_AB*sin(pos(:,1))+r_BC*sin(pos(:,2))-y_c;
    x_c+r_CD*cos(pos(:,2))+r_DE*cos(pos(:,3))-pos(:,4);
    y_c+r_CD*sin(pos(:,2))+r_DE*sin(pos(:,3))-r_FO];
   
    
    end 

%velocity

    function F = myfun1(vel,theta,pos)
r_AB = 5.160;
r_CD = 1.12;
r_BC = 1.12;
r_DE = 2.95;

% Defined
r_EF = 2.5;
r_OA = 0.5;
r_FO = 3;
x_c  = 5;
y_c  = 2;

theta_dot = 50;
F = [-r_OA*theta_dot*sin(theta)-r_AB*vel(1)*sin(pos(:,1))-r_BC*vel(2)*sin(pos(:,2));
    r_OA*theta_dot*cos(theta)+r_AB*vel(1)*cos(pos(:,1))+r_BC*vel(2)*cos(pos(:,2));
    -r_CD*vel(2)*sin(pos(:,2))-r_DE*vel(3)*sin(pos(:,3))-vel(4);
    r_CD*vel(2)*cos(pos(:,2))+r_DE*vel(3)*cos(pos(:,3))];
    
    end 

%%% acceleration

function F = myfun2(acc,theta,pos,vel)
r_AB = 5.160;
r_CD = 1.12;
r_BC = 1.12;
r_DE = 2.95;

% Defined
r_EF = 2.5;
r_OA = 0.5;
r_FO = 3;
x_c  = 5;
y_c  = 2;

theta_dot = 50;
F = [-r_OA*(theta_dot^2)*cos(theta)-r_AB*acc(1)*sin(vel(:,1))-r_AB*(vel(:,1)^2)*cos(vel(:,1))-r_BC*acc(2)*sin(vel(:,2))-r_BC*(vel(:,2)^2)*cos(vel(:,2));
    -r_OA*(theta_dot^2)*sin(theta)+r_AB*acc(1)*cos(vel(:,1))-r_AB*(vel(:,1)^2)*sin(vel(:,1))+r_BC*acc(2)*cos(vel(:,2))-r_BC*(vel(:,2)^2)*cos(vel(:,2));
    -r_CD*acc(2)*sin(vel(:,2))-r_CD*(vel(:,2)^2)*cos(vel(:,2))-r_DE*acc(3)*sin(pos(:,3))-r_DE*(vel(:,3)^2)*cos(pos(:,3))-acc(4);
    r_CD*acc(2)*cos(vel(:,2))-r_CD*(vel(:,2)^2)*sin(vel(:,2))+r_DE*acc(3)*cos(pos(:,3))-r_DE*(vel(:,3)^2)*sin(pos(:,3))];

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
