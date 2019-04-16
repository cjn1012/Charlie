clear all; close all; clc

%Plastic Properties Given 

Mod_El=230000; %Modulus of Elasticity E = 230000psi
Sig_y=2500; %Yeild Strength Sig_y = 2500psi
Sig_uts=3000; %Ultamate Tensile Strength Sig_uts = 3000psi
Density=0.0376; %Density = 0.0376lb/in^3 

%Extra Definitions
w=50; %Omega = 50rpm
theta=(0:1:360).*(pi/180); %Singular Rotation 
k=1; %lbf/in 

%Lengths Determined
r_AB = 4.68;
r_BD = 3.5;
r_DE = 3.4;

%Lengths Defined
r_OA = 0.5;
r_EF = 2.5;

%Geometry of Each Segment 
    %Base = b
    %Height = h 
 
%Base
dimb = .35;
dimh = .30;

b_AB=dimb;
b_DE=dimb;
b_EF=dimb;
%Height

h_AB=dimh;
h_EF=dimh;
h_DE=dimh;

% member BD 3
dim2b = .5;
dim2h = .4;
h_BD=dim2h;
b_BD=dim2b;

%Member OA
dim3b = .35;
dim3h = .15;
b_OA=dim3b;
h_OA=dim3h;


%Dimentions [In Inches]
Width=dimh/2;
Width2=dim2h/2;
Width3=dim3h/2;
Thickness=0.125;
Area=[b_OA*h_OA,b_AB*h_AB,b_BD*h_BD,b_EF*h_EF];

%Y values for Bending Calculations 
y_OA=h_OA/2;
y_AB=h_AB/2;
y_BD=h_BD/2;
y_DE=h_DE/2;


%Inertia Calculations 
I_OA=(b_OA*(h_OA^3))/12;
I_AB=(b_AB*(h_AB^3))/12;
I_BD=(b_BD*(h_BD^3))/12;
I_DE=(b_DE*(h_DE^3))/12;
    %Inertia Matrix 
I=[I_OA,I_AB,I_BD,I_DE];

%Moment Table from Deliverable 1 Solution
    %Values in lb*in
Mmax_OA=0.56;
Mmin_AB=0.00;
Mmax_BD=2.67;
Mmin_DE=0.00;
Mmin_OA=0.54;
Mmax_AB=0.00;
Mmin_BD=0.72;
Mmax_DE=0.00;

Mmin_Bending=[Mmin_OA,Mmin_AB,Mmin_BD,Mmin_DE];
Mmax_Bending=[Mmax_OA,Mmax_AB,Mmax_BD,Mmax_DE];

Min_Bending_Stress=zeros(1,4);
Max_Bending_Stress=zeros(1,4);
for i=1:4
    if i == 3
        Min_Bending_Stress(i)=(Mmin_Bending(i)*(Width2/2)/I(i));
        Max_Bending_Stress(i)=(Mmax_Bending(i)*(Width2/2)/I(i));
    end
    
    Min_Bending_Stress(i)=(Mmin_Bending(i)*(Width/2)/I(i));
    Max_Bending_Stress(i)=(Mmax_Bending(i)*(Width/2)/I(i));
end

%Axial Forces Table from Deliverable 1 Solution
    %Values in lb 
Ax_Mmax_OA=0.4;
Ax_Mmax_AB=1.49;
Ax_Mmin_BD=0.06;
Ax_Mmax_DE=1.50;
Ax_Mmin_OA=0.39;
Ax_Mmin_AB=0.53;
Ax_Mmax_BD=0.36;
Ax_Mmin_DE=0.54;

% %Maxiumum Compression Table from Deliverable 1 Solution
%     %Values in lb
% CMax_OA=-1.49;
% CMax_AB=-1.49;
% CMax_BD=-0.46; %Located in BC Region 
% CMax_DE=-1.5;

%Maxiumum Force in each Pin Table from Deliverable 1 Solutions 
    %Values in lb
FMax_A=1.49;
FMax_B=1.49;
FMax_C=2.97;
FMax_D=1.5;
FMax_E=1.5;
F = [FMax_A,FMax_B,FMax_C,FMax_D,FMax_E]
%Bending Calculations 
    %Bending values calulated by (M*y)/I
    %Maximum
Bending_Max_OA=(Mmax_OA*y_OA)/I_OA;
Bending_Max_AB=(Mmax_AB*y_AB)/I_AB;
Bending_Max_BD=(Mmax_BD*y_BD)/I_BD;
Bending_Max_DE=(Mmax_DE*y_DE)/I_DE;
    %Minimum 
Bending_Min_OA=(Mmin_OA*y_OA)/I_OA;
Bending_Min_AB=(Mmin_AB*y_AB)/I_AB;
Bending_Min_BD=(Mmin_BD*y_BD)/I_BD;
Bending_Min_DE=(Mmin_DE*y_DE)/I_DE;
    %Bending Stress Matrix 
%Max_Bending_Stress=[Bending_Max_OA,Bending_Max_AB,Bending_Max_BD,Bending_Max_DE];
%Min_Bending_Stress=[Bending_Min_OA,Bending_Min_AB,Bending_Min_BD,Bending_Min_DE];

%Axial Calculations 
    %Axial Stress Matrix 
Max_Axial=[Ax_Mmax_OA,Ax_Mmax_AB,Ax_Mmax_BD,Ax_Mmax_DE];
Min_Axial=[Ax_Mmin_OA,Ax_Mmin_AB,Ax_Mmin_BD,Ax_Mmin_DE];

Max_Axial_Stress=Max_Axial./Area;
Min_Axial_Stress=Min_Axial./Area;

%Combined Stress Calculations 
    Max_Comb_Stress=zeros(1,4);
    Min_Comb_Stress=zeros(1,4);
for i=1:4
    Max_Comb_Stress(i)=Max_Bending_Stress(i)+Max_Axial_Stress(i);
    Min_Comb_Stress(i)=Min_Bending_Stress(i)+Min_Axial_Stress(i);
end 

%Mean and Alternating Stress 
Mean_Stress=(Max_Comb_Stress+Min_Comb_Stress)./2;
Alt_Stress=(Max_Comb_Stress-Min_Comb_Stress)./2;

%Buckling & Deflection 
L=[r_OA,r_AB,r_BD,r_DE];
Buckling = (pi^2*Mod_El*I)./(L.^2);

Max_Deflection=L./360;
Actual_Deflection=(Mmax_Bending.*(L.^2))./(9*sqrt(3)*Mod_El.*I);

%Fatigue 
Sf_prime=Sig_uts/2;
Cload=1;
Csize=1;
a_surf=4.511;
b_surf=-0.265;
Csurf=a_surf*Sig_uts^b_surf;
Creliab=0.814;
Ctemp=1;
Sf=Cload*Csize*Csurf*Ctemp*Creliab*Sf_prime;

%Factor of Safety Calculations 
    %Fatigue 
FOS_Fatigue_N1=(Sf*Sig_uts)./((Alt_Stress.*Sig_uts)+(Mean_Stress.*Sf));

    %Yeild
FOS_Yeild_Tension=Sig_y./abs(Alt_Stress - Mean_Stress);
FOS_Yeild_N_2=Sf./Alt_Stress;

FOS_Yeild_N_4=Sig_y./(Alt_Stress + Mean_Stress);

    %Deflection
FOS_Deflection=Max_Deflection./Actual_Deflection;

    %Buckling
FOS_Buckling=Buckling./Max_Axial_Stress;

%Shear Calculations 

%Pin Values and Calculations 
    %Diameter = d [Values in Inches]
d_A=0.137795;
d_B=0.137795;
d_C=0.21;
d_D=0.137795;
d_E=0.137795;
Length=0.19685; 
d = [d_A,d_B,d_C,d_D,d_E]

%Shear Stress Calculations per Pin
Shear_Stress_A=sqrt(3)*FMax_A/(pi*d_A^2)/4;
Shear_Stress_B=sqrt(3)*FMax_B/(pi*d_B^2)/4;
Shear_Stress_C=sqrt(3)*FMax_C/(pi*d_C^2)/4;
Shear_Stress_D=sqrt(3)*FMax_D/(pi*d_D^2)/4;
Shear_Stress_E=sqrt(3)*FMax_E/(pi*d_E^2)/4;
%Shear Stress Array
Shear_Stress=[Shear_Stress_A,Shear_Stress_B,Shear_Stress_C,Shear_Stress_D,Shear_Stress_E];

%Bearing Stress Calculations per Pin
Bearing_Stress_A=FMax_A/((pi*d_A*Length)/4);
Bearing_Stress_B=FMax_B/((pi*d_B*Length)/4);
Bearing_Stress_C=FMax_C/((pi*d_C*Length)/4);
Bearing_Stress_D=FMax_D/((pi*d_D*Length)/4);
Bearing_Stress_E=FMax_E/((pi*d_E*Length)/4);
%Bearing Stress Array
Bearing_Stress=[Bearing_Stress_A,Bearing_Stress_B,Bearing_Stress_C,Bearing_Stress_D,Bearing_Stress_E];

%Factor of Safety Pin Calculations
FOS_Shear=Sig_y./Shear_Stress;
FOS_Bearing=Sig_y./Bearing_Stress;

%Tearout Calculations 
Tearout_Stress_OA=sqrt(3)*Max_Axial(1)/(2*b_OA*h_OA);
Tearout_Stress_AB=sqrt(3)*Max_Axial(2)/(2*b_AB*h_AB);
Tearout_Stress_BD=sqrt(3)*Max_Axial(3)/(2*b_BD*h_BD);
Tearout_Stress_DE=sqrt(3)*Max_Axial(4)/(2*b_DE*h_DE);
%Tearout Array
Tearout_Stress=[Tearout_Stress_OA,Tearout_Stress_AB,Tearout_Stress_BD,Tearout_Stress_DE];

%Factor of Safety for Tearout
FOS_Tearout=Sig_y./Tearout_Stress;



% Creating FOS Table
t = table(FOS_Fatigue_N1',FOS_Yeild_Tension',FOS_Yeild_N_2'...
    ,FOS_Yeild_N_4',FOS_Deflection',FOS_Buckling',FOS_Tearout');

t2 = table(FOS_Bearing',FOS_Shear');







