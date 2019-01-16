%THIS CODE IS FOR GENERATING THE blockMeshDict FILE FOR A TUNNEL SIMULATION
%INPUT (yplus - yplus size, er - expansion ratio, xs and ys cells in
%streamwise direction (for contraction) and azimuthal direction (for 45
%degrees), gpx and rpx - number of spline points
%REQUIRES Rratio.m
%NOTE: GEOMETRY MAY BE SLIGHTLY DIFFERENT SINCE BLOCKMESH IS USED RATHER
%THAN THE ORIGINAL STL FILE; INLET IS 1/4 D UPSTREAM, OUTLET 10 INCHES
%DOWNSTREAM OF REAL TEST SECTION

%SCALE
a=1;

%GEOMETRY POINTS IN x and r
gpx=1000;
rpx=200;

%Y plus, expansion ratios, etc.
yplus=1e-3;     %value of y+=1
er=1.1;         %expansion ratio
R=Rratio(3*a,er,yplus);
nz=ceil(1+log(R)/log(er));

%GRADING AND CELLS
xs=30;  %along length of contraction %100
ys=3;  %20
zs=nz;  %50 %or calculated above
xs2=floor(xs*36/27);    %along length of ts
xg=1;
yg=1;   
zg=R;   %10 %or calulated above

%POINTS
p0=[0 0 0];
p1=[27*a 0 0];
p2=[27*a 3*a 0];
p3=[0 9*a 0];
p4=[0 sqrt(2)/2*(9*a) sqrt(2)/2*(9*a)];
p5=[27*a 3*a 3*a];
p6=[27*a 0 3*a];
p7=[0 0 9*a];
p8=[(27+36)*a 0 0];
p9=[(27+36)*a 0 3*a];
p10=[(27+36)*a 3*a 0];
p11=[(27+36)*a 3*a 3*a];

L=p2(1)-p0(1);

yin=9*a;
yout=2/sqrt(2)*3*a;
for i=1:gpx
    x(i)=L/(gpx+1)*i;
    y(i)=sqrt(2)/2*(yin-(yin-yout)*(6*(x(i)/L)^5-15*(x(i)/L)^4+10*(x(i)/L)^3));
end
z=y;
    
yin=9*a;
yout=3*a;
for i=1:gpx
    x2(i)=L/(gpx+1)*i;
    y2(i)=yin-(yin-yout)*(6*(x2(i)/L)^5-15*(x2(i)/L)^4+10*(x2(i)/L)^3);
    z2(i)=0;
end

for i=1:rpx
    x3(i)=0;
    xtmp=pi()/4/(rpx+1)*i;
    y3(i)=p3(2)*cos(xtmp);
    z3(i)=p3(2)*sin(xtmp);
end

for i=1:rpx
    x4(i)=0;
    xtmp=pi()/4+pi()/4/(rpx+1)*i;
    y4(i)=p3(2)*cos(xtmp);
    z4(i)=p3(2)*sin(xtmp);
end

yin=9*a;
yout=3*a;
for i=1:gpx
    x5(i)=L/(gpx+1)*i;
    y5(i)=0;
    z5(i)=yin-(yin-yout)*(6*(x5(i)/L)^5-15*(x5(i)/L)^4+10*(x5(i)/L)^3);    
end



figure(1)
plot3(x,y,z)
hold on
plot3(x2,y2,z2)
plot3(x3,y3,z3)
plot3(x4,y4,z4)
plot3(x5,y5,z5)
title('spline profiles')


f=fopen('blockMeshDict','w');

    %WRITING PREAMBLE
    fprintf(f, '/*--------------------------------*- C++ -*----------------------------------*\\ \n');
    fprintf(f, '|  =========                 |                                                | \n');
    fprintf(f, '|  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox          | \n');
    fprintf(f, '|   \\\\    /   O peration     | Version:  2.2.1                                | \n');
    fprintf(f, '|    \\\\  /    A nd           | Web:      www.OpenFOAM.org                     | \n');
    fprintf(f, '|     \\\\/     M anipulation  |                                                | \n');
    fprintf(f, '\\*---------------------------------------------------------------------------*/ \n');
    fprintf(f, '\n');
    fprintf(f, 'FoamFile \n');
    fprintf(f, '{ \n');
    fprintf(f, '    version     2.0; \n');
    fprintf(f, '    format      ascii; \n');
    fprintf(f, '    class       dictionary; \n');
    fprintf(f, '    object      blockMeshDict; \n');
    fprintf(f, '} \n');
    fprintf(f, ' \n');
    fprintf(f, 'convertToMeters 0.0254; \n');
    fprintf(f, ' \n');

    %VERTICES
    fprintf(f, 'vertices \n');
    fprintf(f, '( \n'); 

    fprintf(f, '(%f %f %f) //point 0 \n', p0(1), p0(2), p0(3)); 
    fprintf(f, '(%f %f %f) //point 1 \n', p1(1), p1(2), p1(3)); 
    fprintf(f, '(%f %f %f) //point 2 \n', p2(1), p2(2), p2(3)); 
    fprintf(f, '(%f %f %f) //point 3 \n', p3(1), p3(2), p3(3)); 
    fprintf(f, '(%f %f %f) //point 4 \n', p4(1), p4(2), p4(3)); 
    fprintf(f, '(%f %f %f) //point 5 \n', p5(1), p5(2), p5(3)); 
    fprintf(f, '(%f %f %f) //point 6 \n', p6(1), p6(2), p6(3)); 
    fprintf(f, '(%f %f %f) //point 7 \n', p7(1), p7(2), p7(3)); 
    fprintf(f, '(%f %f %f) //point 8 \n', p8(1), p8(2), p8(3)); 
    fprintf(f, '(%f %f %f) //point 9 \n', p9(1), p9(2), p9(3)); 
    fprintf(f, '(%f %f %f) //point 10 \n', p10(1), p10(2), p10(3)); 
    fprintf(f, '(%f %f %f) //point 11 \n', p11(1), p11(2), p11(3)); 
    
    fprintf(f, '); \n'); 
    
    %BLOCKS
    fprintf(f, 'blocks \n');
    fprintf(f, '( \n'); 
    
    fprintf(f, 'hex (3 2 5 4 0 1 1 0) (%u %u %u) simpleGrading (%f %f %f) \n',xs, ys, zs, xg, yg, zg);
    fprintf(f, 'hex (4 5 6 7 0 1 1 0) (%u %u %u) simpleGrading (%f %f %f) \n',xs, ys, zs, xg, yg, zg);
    fprintf(f, 'hex (2 10 11 5 1 8 8 1) (%u %u %u) simpleGrading (%f %f %f) \n',xs2, ys, zs, xg, yg, zg);
    fprintf(f, 'hex (5 11 9 6 1 8 8 1) (%u %u %u) simpleGrading (%f %f %f) \n',xs2, ys, zs, xg, yg, zg);
    
    fprintf(f, '); \n');
    
    %EDGES
    fprintf(f, 'edges \n');
    fprintf(f, '( \n');

    fprintf(f, 'spline 4 5 \n');
    fprintf(f, '( \n');
    for j=1:length(x)
    fprintf(f, '(%f %f %f) \n', x(j), y(j), z(j));
    end
    fprintf(f, ') \n');
        
    fprintf(f, 'spline 3 2 \n');
    fprintf(f, '( \n');
    for j=1:length(x2)
        fprintf(f, '(%f %f %f) \n', x2(j), y2(j), z2(j));
    end
    fprintf(f, ') \n');
    
    fprintf(f, 'spline 3 4 \n');
    fprintf(f, '( \n');
    for j=1:length(x3)
        fprintf(f, '(%f %f %f) \n', x3(j), y3(j), z3(j));
    end
    fprintf(f, ') \n');
    
    fprintf(f, 'spline 4 7 \n');
    fprintf(f, '( \n');
    for j=1:length(x4)
        fprintf(f, '(%f %f %f) \n', x4(j), y4(j), z4(j));
    end
    fprintf(f, ') \n');
    
    fprintf(f, 'spline 7 6 \n');
    fprintf(f, '( \n');
    for j=1:length(x5)
        fprintf(f, '(%f %f %f) \n', x5(j), y5(j), z5(j));
    end
    fprintf(f, ') \n');
    
    fprintf(f, '); \n');    
    
    %BOUNDARIES
    fprintf(f, 'boundary \n');
    fprintf(f, '( \n');  

    %Inlet
    fprintf(f, '  Inlet \n');
    fprintf(f, '    { \n');
    fprintf(f, '      type patch; \n');
    fprintf(f, '      faces \n');
    fprintf(f, '      ( \n');   
    
    fprintf(f, '       (4 0 0 7) \n');   
    fprintf(f, '       (3 0 0 4) \n'); 
    
    fprintf(f, '      ); \n'); 
    fprintf(f, '    } \n');  
    
    %Outlet
    fprintf(f, '  Outlet \n');
    fprintf(f, '    { \n');
    fprintf(f, '      type patch; \n');
    fprintf(f, '      faces \n');
    fprintf(f, '      ( \n');   
    
    fprintf(f, '       (11 9 8 8) \n');   
    fprintf(f, '       (10 11 8 8) \n'); 
    
    fprintf(f, '      ); \n'); 
    fprintf(f, '    } \n');  
    
	%Bottom
    fprintf(f, '  Bottom \n');
    fprintf(f, '    { \n');
    fprintf(f, '      type symmetryPlane; \n');
    fprintf(f, '      faces \n');
    fprintf(f, '      ( \n');   
    
    fprintf(f, '       (3 2 1 0) \n');   
    fprintf(f, '       (2 10 8 1) \n'); 
    
    fprintf(f, '      ); \n'); 
    fprintf(f, '    } \n');  
    
    %Side
    fprintf(f, '  Side \n');
    fprintf(f, '    { \n');
    fprintf(f, '      type symmetryPlane; \n');
    fprintf(f, '      faces \n');
    fprintf(f, '      ( \n');   
    
    fprintf(f, '       (0 1 6 7) \n');   
    fprintf(f, '       (1 8 9 6) \n'); 
    
    fprintf(f, '      ); \n'); 
    fprintf(f, '    } \n');  
    
	%Walls
    fprintf(f, '  Walls \n');
    fprintf(f, '    { \n');
    fprintf(f, '      type wall; \n');
    fprintf(f, '      faces \n');
    fprintf(f, '      ( \n');   
    
    fprintf(f, '       (4 5 2 3) \n');   
    fprintf(f, '       (7 6 5 4) \n'); 
    fprintf(f, '       (5 11 10 2) \n');   
    fprintf(f, '       (9 11 5 6) \n'); 
    
    fprintf(f, '      ); \n'); 
    fprintf(f, '    } \n');  
    
% 	%Empties    NOT NEEDED BECAUSE OF FACE MATCHING! :)
%     fprintf(f, '  Empties \n');
%     fprintf(f, '    { \n');
%     fprintf(f, '      type empty; \n');
%     fprintf(f, '      faces \n');
%     fprintf(f, '      ( \n');   
%     
%     fprintf(f, '       (0 1 1 0) \n');   
%     
%     fprintf(f, '      ); \n'); 
%     fprintf(f, '    } \n');  
    
    fprintf(f, '); \n'); 
    
    fprintf(f, 'mergePatchPairs \n'); 
    fprintf(f, '( \n'); 
    fprintf(f, '); \n'); 
    fclose(f);
    
%FILE CLOSED

numberOfCells=(xs+xs2)*ys*zs*2








