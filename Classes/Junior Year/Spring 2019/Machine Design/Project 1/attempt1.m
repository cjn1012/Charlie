clear all; close all;

theta = (0:1:360).*(pi/180);
for i= 1:length(theta)
    x0 = [5.7*pi/180,120*pi/180,0*pi/180,2.5];        
    options=optimset('Display','iter');  
    x(i,:) = fsolve(@(x) myfun(x,theta(i)),x0,options); 
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
    x_c+r_CD*sin(x(:,2))+r_DE*cos(pi+x(:,3))-r_EF;
    y_c+r_CD*cos(x(:,2))+r_DE*sin(pi+x(:,3))-x(:,4)];
   
    
    end %THIS END ISNT WORKING FIX ME
