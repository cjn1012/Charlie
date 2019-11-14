clear all
close all

LVT_Acc = xlsread('LVTandAccel.csv', 'A8:A100007');
Acc = xlsread('LVTandAccel.csv', 'B8:B100007');
Time = linspace(-.34745,-.34745+0.00001168*100000,100000)';
%% a
figure(1)
plot(Time,LVT_Acc)
hold on
plot(Time,Acc)
xlabel('Time (s)')
ylabel('Voltage (V)')
legend('LVT Output','Accelerometer Output')


Acc_Fall = -.09; %Volts for 1g
Sens_Acc = abs(Acc_Fall/32.2) %volts/in/sec2

% At t = -.06, object starts falling and is constant till t=-.04
for x = 1:100000
    if Time(x)>-.06
        t1_x = x;
        break
    end
end

for x = 1:100000
    if Time(x)>-.04
        t2_x = x;
        break
    end
end

Change_V = 32.2*.02;
Delta_Voltage = abs(LVT_Acc(t2_x) - LVT_Acc(t1_x));
Sens_LVT_Acc = Delta_Voltage/Change_V %volt/in/sec
%%
%b

Acc_In = Acc./Sens_Acc;
Acc_In = Acc_In - mean(Acc_In(1:1500));
LVT_Acc_In = LVT_Acc./Sens_LVT_Acc;


Int_Acc = cumtrapz(Time,Acc_In);
figure(2)
plot(Time,Int_Acc,Time,LVT_Acc_In)
xlabel('Time (s)')
ylabel('Velocity (in/s)')


%% c

LVT_Force = xlsread('LVTandForce.csv', 'A8:A100007');
Force = xlsread('LVTandForce.csv', 'B8:B100007');
Force = Force-mean(Force(1:100));
%beginning of the drop
for x = 1:100000
    if Acc_In(x)<-5 
        x_drop = x;
        break
    end
end

%first contact of core and foam
for x = 1:100000
    if Force(x)> .03
        x_hit = x;
        break
    end
end

%max displacement --- max force!

Force_Max = max(Force);
for x = 1:100000
    if Force(x)== Force_Max
        x_maxdisp = x;
        break
    end
end

% finding first bounce
th = 0.03;
[peakLoc, peakMag] = peakfinder(Force, th,'minima',-1);
figure(3)
plot(Time,Force,Time(peakLoc),Force(peakLoc))

x_firstbounce = peakLoc(2); %possible, might not be a bounce, very close

%permanent contact
% once it hits foam, it doesnt seem to bounce at all


figure(4)
plot(Time,LVT_Acc)
hold on
plot(Time,Acc)
plot(Time(x_drop),LVT_Acc(x_drop),'o',Time(x_hit),LVT_Acc(x_hit),'o',Time(x_maxdisp),LVT_Acc(x_maxdisp),'o',Time(x_firstbounce),LVT_Acc(x_firstbounce),'o')
text(.25,-.5,strcat('Beg. of Drop is ' ,{' '},num2str(Time(x_drop)),{' '}, 'sec'))
text(.25,-.75,strcat('First mass contact is' ,{' '},num2str(Time(x_hit)),{' '}, 'sec'))
text(.25,-1,strcat('Max disp. is' ,{' '},num2str(Time(x_maxdisp)),{' '}, 'sec'))
text(.25,-1.25,strcat('Possible first bounce is ' ,{' '},num2str(Time(x_firstbounce)),{' '}, 'sec'))
text(.25,-1.5,'Permanent contact after first bounce')

xlabel('Time (s)')
ylabel('Voltage (V)')
legend('LVT Output','Accelerometer Output')

%% d. maximum velocity of core
LVT_Acc_In;

Max_Vel = max(abs(LVT_Acc_In)) %in/s

%%
figure(5)
plot(Time,LVT_Force,Time,Force)
xlabel('Time (s)')
ylabel('Voltage (V)')



Force_lbf = Force./.491; %lbf

for x = 1:100000
    if abs(LVT_Acc_In(x))== Max_Vel
        x_maxvel = x;
        break
    end
end

Force_maxvel = Force_lbf(x_maxvel)
Force_steady = mean(Force_lbf(end-1000,end))

total_mass = Force_steady/32.2















