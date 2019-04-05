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

theta_dot = 50*2*pi/60;


%Loop #1
R_OAx=(r_OA/2)*cos(theta);
R_OAy=(r_OA/2)*sin(theta);
R_ABx=(r_AB/2)*cos(pos(:,1));
R_ABy=(r_AB/2)*sin(pos(:,1));
R_BCx=(r_BC/2)*cos(pos(:,2));
R_BCy=(r_BC/2)*sin(pos(:,2));


%Loop #2
R_CDx=(r_CD/2)*sin(pos(:,2));
R_CDy=(r_CD/2)*cos(pos(:,2));
R_DEx=(r_DE/2)*cos(pos(:,3));
R_DEy=(r_DE/2)*sin(pos(:,3));
R_EFx=pos(:,4);
R_EFy=3;
R_FOx=0;
R_FOy=1.5;


% Force Calculations

Density = .0376; %lb/in^3
Thickness = 0.25;
Width = 1.0;

m_OA = Density*Thickness*Width*r_OA;
m_AB = Density*Thickness*Width*r_AB;
m_BC = Density*Thickness*Width*r_BC;
m_CD = Density*Thickness*Width*r_CD;
m_DE = Density*Thickness*Width*r_DE;


Unknowns = 14;
for i = 1:length(theta)
    x0 = zeros(Unknowns,1)';
    options=optimset('Display','iter');
    w(i,:) = fsolve(@(w) myfun3(w,pos(i,2),pos(i,4),acc(:,4),...
        R_OAx(i),R_OAy(i),R_ABx(i),R_ABy(i),R_BCx(i),R_BCy(i),R_CDx(i),R_CDy(i))...
        ,x0,options);

end



%% Graphs

r_EFy = r_FO;
yEF(1:361) = r_EFy;
figure(1);
plot(pos(:,4),yEF, 'Linewidth', 1)
hold on
plot(.5*cos(theta),.5*sin(theta),'Linewidth', 1)
plot(1.12*cos(pos(:,2))+x_c , 1.12*sin(pos(:,2))+y_c, 'Linewidth', 1)
plot(1.12*cos(pos(:,2)+pi)+x_c , 1.12*sin(pos(:,2)+pi)+y_c, 'Linewidth', 1)
grid minor;
xlabel(' X Position [in] ')
ylabel(' Y Position [in] ')
title('Trajectories of Points')
legend('Point E','Point A','Point D','Point B','Location','northwest')
axis([-1 6 -1 3.5])

figure(2);
title('X Position of Points A and E Vs. Crank Angle')
hold on
plot(theta*(180/pi),.5*cos(theta),theta*(180/pi),pos(:,4))
xlabel(' Crank Angle [°]  ')
ylabel(' X Position [in] ')
title('Trajectories of Points')
legend('Point A','Point E','Location','northwest')
xlim([0 360])
hold off

figure(3)
hold on
title('X-Component Linear Velocities of Points A and E Vs. Crank Angle')
plot(theta*(180/pi),-.5*sin(theta)*theta_dot)
plot(theta*(180/pi),vel(:,4))
xlabel(' Crank Angle [°]  ')
ylabel(' \nu_x [in/s] ')
legend('Point A','Point E','Location','northwest')
xlim([0 360])
hold off

figure(4)
hold on
title('X-Component Linear Accerlation Vs. Crank Angle')
plot(theta*(180/pi),-.5*cos(theta))
plot(theta*(180/pi),acc(:,4))
plot(theta*(180/pi),pos(:,1)*r_OA/2)
plot(theta*(180/pi),pos(:,2)*r_AB/2)
plot(theta*(180/pi),pos(:,3)*r_BC/2)
plot(theta*(180/pi),pos(:,4)*r_CD/2)
xlabel(' Crank Angle [°]  ')
ylabel(' \alpha_x [in/s^2] ')
legend('Point A','Point E','Member 1','Member 2','Member 3','Member 4','position','southeast')
xlim([0 360])
hold off

for i = 1:length(theta)
    F_Joint_O(i) = pythag(w(i,1),w(i,3));
    F_Joint_A(i) = pythag(w(i,2),w(i,4));
    F_Joint_B(i) = pythag(w(i,6),w(i,7));
    F_Joint_C(i) = pythag(w(i,9),w(i,8));
    F_Joint_D(i) = pythag(w(i,12),w(i,13));
    F_Joint_E(i) = pythag(w(13)*cos(pos(i,3)),(2-pos(i,4)));
    
end

figure(5)
hold on
title('Forces in Joints')
plot(theta*(180/pi),F_Joint_O)
plot(theta*(180/pi),F_Joint_A)
plot(theta*(180/pi),F_Joint_B)
plot(theta*(180/pi),F_Joint_C)
plot(theta*(180/pi),F_Joint_D)
plot(theta*(180/pi),F_Joint_E)
xlabel(' Crank Angle [°]  ')
ylabel(' Magnitude of Force [lb] ')
xlim([0 360])
legend('Joint O','Joint A','Joint B','Joint C','Joint D','Joint E')

%%  Axial and Shear Force

length_member_OA = linspace(0,r_OA,361);
length_member_AB = linspace(0,r_AB,361);
length_member_BD = linspace(0,2*r_BC,361);
length_member_DE = linspace(0,r_DE,361);

for i = 1:length(theta)
    for j = 1:length(length_member_OA)
        axial_OA(i,j) = BigBootyBitchz(w(i,1),w(i,3),w(i,2),w(i,4));
        shear_OA(i,j) = real(BigBootyBitchz2(w(i,1),w(i,3),w(i,2),w(i,4)));

    end
end


for i = 1:length(theta)
    for j = 1:length(length_member_AB)
        axial_AB(i,j) = BigBootyBitchz(w(i,2),w(i,4),w(i,6),w(i,7));
        shear_AB(i,j) = real(BigBootyBitchz2(w(i,2),w(i,4),w(i,6),w(i,7)));
    end
end

for i = 1:length(theta)
    for j = 1:length(length_member_BD)
        axial_BD(i,j) = BigBootyBitchz(w(i,6),w(i,7),w(i,12),w(i,13));
        shear_BD(i,j) = real(BigBootyBitchz2(w(i,6),w(i,7),w(i,12),w(i,13)));
    end
end

for i = 1:length(theta)
    for j = 1:length(length_member_DE)
        axial_DE(i,j) = BigBootyBitchz(w(i,12),w(i,13),w(13)*cos(pos(i,3)),(2-pos(i,4)));
        shear_DE(i,j) = real(BigBootyBitchz2(w(i,12),w(i,13),w(13)*cos(pos(i,3)),(2-pos(i,4))));
    end
end




[Max_Axial_OA,Max_Axial_OA_L] = max(max(axial_OA));
[Max_Axial_AB,Max_Axial_AB_L] = max(max(axial_AB));
[Max_Axial_BD,Max_Axial_BD_L] = max(max(axial_BD));
[Max_Axial_DE,Max_Axial_DE_L] = max(max(axial_DE));
[Max_Shear_OA,Max_Shear_OA_L] = max(max(shear_OA));
[Max_Shear_AB,Max_Shear_AB_L] = max(max(shear_AB));
[Max_Shear_BD,Max_Shear_BD_L] = max(max(shear_BD));
[Max_Shear_DE,Max_Shear_DE_L] = max(max(shear_DE));

[Min_Axial_OA,Min_Axial_OA_L] = min(min(axial_OA));
[Min_Axial_AB,Min_Axial_AB_L] = min(min(axial_AB));
[Min_Axial_BD,Min_Axial_BD_L] = min(min(axial_BD));
[Min_Axial_DE,Min_Axial_DE_L] = min(min(axial_DE));
[Min_Shear_OA,Min_Shear_OA_L] = min(min(shear_OA));
[Min_Shear_AB,Min_Shear_AB_L] = min(min(shear_AB));
[Min_Shear_BD,Min_Shear_BD_L] = min(min(shear_BD));
[Min_Shear_DE,Min_Shear_DE_L] = min(min(shear_DE));





% Axial

figure(6);
surf(length_member_OA, theta*180/pi, axial_OA,'edgecolor' , 'none', 'facecolor','interp')
colorbar
xlabel('Length [in]')
ylabel('Crank Angle [°]')
zlabel('Axial Force [lb]')
title('Axial Force on Member OA')
ylim([0 360])

figure(7);
surf(length_member_AB, theta*180/pi, axial_AB,'edgecolor' , 'none', 'facecolor','interp')
colorbar
xlabel('Length [in]')
ylabel('Crank Angle [°]')
zlabel('Axial Force [lb]')
title('Axial Force on Member AB')
ylim([0 360])

figure(8);
surf(length_member_BD, theta*180/pi, axial_BD,'edgecolor' , 'none', 'facecolor','interp')
colorbar
xlabel('Length [in]')
ylabel('Crank Angle [°]')
zlabel('Axial Force [lb]')
title('Axial Force on Member BD')
ylim([0 360])

figure(9);
surf(length_member_DE, theta*180/pi, axial_DE,'edgecolor' , 'none', 'facecolor','interp')
colorbar
xlabel('Length [in]')
ylabel('Crank Angle [°]')
zlabel('Axial Force [lb]')
title('Axial Force on Member DE')
ylim([0 360])



figure(10);
surf(length_member_OA, theta*180/pi, shear_OA,'edgecolor' , 'none', 'facecolor','interp')
colorbar
xlabel('Length [in]')
ylabel('Crank Angle [°]')
zlabel('Shear Force [lb]')
title('Shear Force on Member OA')
ylim([0 360])


figure(11);
surf(length_member_AB, theta*180/pi, shear_AB,'edgecolor' , 'none', 'facecolor','interp')
colorbar
xlabel('Length [in]')
ylabel('Crank Angle [°]')
zlabel('Shear Force [lb]')
title('Shear Force on Member AB')
ylim([0 360])

figure(12);
surf(length_member_BD, theta*180/pi, shear_BD,'edgecolor' , 'none', 'facecolor','interp')
colorbar
xlabel('Length [in]')
ylabel('Crank Angle [°]')
zlabel('Shear Force [lb]')
title('Shear Force on Member BD')
ylim([0 360])


figure(13);
surf(length_member_DE, theta*180/pi, shear_DE,'edgecolor' , 'none', 'facecolor','interp')
colorbar
xlabel('Length [in]')
ylabel('Crank Angle [°]')
zlabel('Shear Force [lb]')
title('Shear Force on Member DE')
ylim([0 360])




%% Moment



Moment_OA = shear_OA.*length_member_OA;
Moment_AB = shear_AB.*length_member_AB;
Moment_BD = shear_BD.*length_member_BD;
Moment_DE = shear_DE.*length_member_DE;


[Max_Moment_OA,Max_Moment_OA_L] = max(Moment_OA);
Max_Moment_OA_L = Max_Moment_OA_L(2);
[Max_Moment_AB,Max_Moment_AB_L] = max(Moment_AB);
Max_Moment_AB_L = Max_Moment_AB_L(2);
[Max_Moment_BD,Max_Moment_BD_L] = max(Moment_BD);
Max_Moment_BD_L = Max_Moment_BD_L(2);
[Max_Moment_DE,Max_Moment_DE_L] = max(Moment_DE);
Max_Moment_DE_L = Max_Moment_DE_L(2);
[Max_Moment_OA2,Max_Moment_OA_L2] = max(Max_Moment_OA);
[Max_Moment_AB2,Max_Moment_AB_L2] = max(Max_Moment_AB);
[Max_Moment_BD2,Max_Moment_BD_L2] = max(Max_Moment_BD);
[Max_Moment_DE2,Max_Moment_DE_L2] = max(Max_Moment_DE);

[Min_Moment_OA,Min_Moment_OA_L] = min(Moment_OA);
Min_Moment_OA_L = Min_Moment_OA_L(2);
[Min_Moment_AB,Min_Moment_AB_L] = min(Moment_AB);
Min_Moment_AB_L = Min_Moment_AB_L(2);
[Min_Moment_BD,Min_Moment_BD_L] = min(Moment_BD);
Min_Moment_BD_L = Min_Moment_BD_L(2);
[Min_Moment_DE,Min_Moment_DE_L] = min(Moment_DE);
Min_Moment_DE_L = Min_Moment_DE_L(2);
[Min_Moment_OA2,Min_Moment_OA_L2] = min(Min_Moment_OA);
[Min_Moment_AB2,Min_Moment_AB_L2] = min(Min_Moment_AB);
[Min_Moment_BD2,Min_Moment_BD_L2] = min(Min_Moment_BD);
[Min_Moment_DE2,Min_Moment_DE_L2] = min(Min_Moment_DE);




figure(14);
surf(length_member_OA, theta*180/pi, Moment_OA,'edgecolor' , 'none', 'facecolor','interp')
colorbar
xlabel('Length [in]')
ylabel('Crank Angle [°]')
zlabel('Internal Bending Moment [lb-in]')
title('Internal Bending Moment on Member OA')
ylim([0 360])

figure(15);
surf(length_member_AB, theta*180/pi, Moment_AB,'edgecolor' , 'none', 'facecolor','interp')
colorbar
xlabel('Length [in]')
ylabel('Crank Angle [°]')
zlabel('Internal Bending Moment [lb-in]')
title('Internal Bending Moment on Member AB')
ylim([0 360])

figure(16);
surf(length_member_BD, theta*180/pi, Moment_BD,'edgecolor' , 'none', 'facecolor','interp')
colorbar
xlabel('Length [in]')
ylabel('Crank Angle [°]')
zlabel('Internal Bending Moment [lb-in]')
title('Internal Bending Moment on Member BD')
ylim([0 360])

figure(17);
surf(length_member_DE, theta*180/pi, Moment_DE,'edgecolor' , 'none', 'facecolor','interp')
colorbar
xlabel('Length [in]')
ylabel('Crank Angle [°]')
zlabel('Internal Bending Moment [lb-in]')
title('Internal Bending Moment on Member DE')
ylim([0 360])

    


Table = [Max_Axial_OA,Max_Axial_AB,Max_Axial_BD,Max_Axial_DE,Max_Shear_OA,Max_Shear_AB,Max_Shear_BD,Max_Shear_DE,...
    Min_Axial_OA,Min_Axial_AB,Min_Axial_BD,Min_Axial_DE,Min_Shear_OA,Min_Shear_AB,Min_Shear_BD,Min_Shear_DE...
    ,Max_Moment_OA2,Max_Moment_AB2,Max_Moment_BD2,Max_Moment_DE2,Min_Moment_OA2,Min_Moment_AB2,Min_Moment_BD2,Min_Moment_DE2]

Table2 = [max(F_Joint_O),max(F_Joint_A),max(F_Joint_B),max(F_Joint_C),max(F_Joint_D),max(F_Joint_E)]'


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

theta_dot = 50*2*pi/60;

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

%     F = [r_OA*cos(theta)+r_AB*cos(pos(:,1))+r_BC*cos(pos(:,2))-x_c;
%     r_OA*sin(theta)+r_AB*sin(pos(:,1))+r_BC*sin(pos(:,2))-y_c;
%     x_c+r_CD*cos(pos(:,2))+r_DE*cos(pos(:,3))-pos(:,4);
%     y_c+r_CD*sin(pos(:,2))+r_DE*sin(pos(:,3))-r_FO];

theta_dot = 50*2*pi/60;
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

% F = [-r_OA*theta_dot*sin(theta)-r_AB*vel(1)*sin(pos(:,1))-r_BC*vel(2)*sin(pos(:,2));
%     r_OA*theta_dot*cos(theta)+r_AB*vel(1)*cos(pos(:,1))+r_BC*vel(2)*cos(pos(:,2));
%     -r_CD*vel(2)*sin(pos(:,2))-r_DE*vel(3)*sin(pos(:,3))-vel(4);
%     r_CD*vel(2)*cos(pos(:,2))+r_DE*vel(3)*cos(pos(:,3))];

theta_dot = 50*2*pi/60;
F = [-r_OA*(theta_dot^2)*cos(theta)-r_AB*acc(1)*sin(pos(:,1))-r_AB*(vel(:,1)^2)*cos(pos(:,1))-r_BC*acc(2)*sin(pos(:,2))-r_BC*(vel(:,2)^2)*cos(pos(:,2));
    -r_OA*(theta_dot^2)*sin(theta)+r_AB*acc(1)*cos(pos(:,1))-r_AB*(vel(:,1)^2)*sin(pos(:,1))+r_BC*acc(2)*cos(pos(:,2))-r_BC*(vel(:,2)^2)*sin(pos(:,2));
    -r_CD*acc(2)*sin(pos(:,2))-r_CD*(vel(:,2)^2)*cos(pos(:,2))-r_DE*acc(3)*sin(pos(:,3))-r_DE*(vel(:,3)^2)*cos(pos(:,3))-acc(4);
    r_CD*acc(2)*cos(pos(:,2))-r_CD*(vel(:,2)^2)*sin(pos(:,2))+r_DE*acc(3)*cos(pos(:,3))-r_DE*(vel(:,3)^2)*sin(pos(:,3))];

end

function F = myfun3(w,x_2,x_4,a_4,...
    R_OAx,R_OAy,R_ABx,R_ABy,R_BCx,R_BCy,R_CDx,R_CDy)


F = [
    %Member 1
    w(1)+w(2);
    w(3)+w(4);
    w(5)+R_OAy*w(4)-R_OAx*w(2)-R_OAy*w(1)+R_OAx*w(3);
    %Member 2
    -w(2)+w(6)*cos(x_2);
    -w(4)+w(6)*sin(x_2);
    -w(2)*R_ABx+w(4)*R_ABy-w(6)*R_ABx+w(7)*R_ABy;
    %Member 3
    w(8)-w(6)+w(9);
    w(10)+w(7)+w(11);
    w(6)*R_BCx+w(7)*R_BCy+w(8)*R_BCx+w(10)*R_BCy;
    %Member 4
    -w(8)+w(12);
    -w(10)+w(13);
    -w(10)*R_CDx+w(8)*R_CDy-w(13)*R_CDx+w(12)*R_CDy;
    %Member E
    (2-x_4)-w(12)-.05*a_4;
    w(14)+w(13);
    ]

end
    
    
function TotalF = pythag(x,y)

TotalF = sqrt(x^2+y^2);

end

    
function axialf = BigBootyBitchz(x1,y1,x2,y2)
angle1 = atan(y1/x1);
angle2 = atan(y2/x2);
axialf = sqrt(pythag(x1,y1)*cos(angle1)+pythag(x2,y2)*cos(angle2));
end
function axialf = BigBootyBitchz2(x1,y1,x2,y2)
angle1 = atan(y1/x1);
angle2 = atan(y2/x2);
axialf = sqrt(pythag(x1,y1)*sin(angle1)+pythag(x2,y2)*sin(angle2));
end

