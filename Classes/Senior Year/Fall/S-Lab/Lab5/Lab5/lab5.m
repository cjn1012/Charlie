%% Part 2
clear all
close all
Input_ei=[1.5 2 3 4 5 6 7 8 9 10]; %Volts
Tach_eo=[.265 .561 1.14 1.72 2.38 2.994 3.579 4.195 4.794 5.405]; %Volts
MUT_em=[.390 .867 1.756 2.661 3.85 4.85 5.806 6.784 7.768 8.751]; %Volts

%a
Ktach=3; %Volts/KRPM

%b
for i=1:length(Tach_eo)
Omega(i)=Tach_eo(i)/Ktach;
end
MUTfit=polyfit(Omega,MUT_em,1);
Ke=MUTfit(1); % Ke motor voltage constant (Volts/KRPM)
Keactual=(4.39+5.37)/2;
figure(1)
bestfit=(MUTfit(1)*Omega)+MUTfit(2);
plot(Omega,MUT_em,'bo-',Omega,bestfit,'r-')
ylabel('MUT Voltage, e_m (Volts)')
xlabel('Angular Velocity, \omega (KRPM)')
title('MUT Voltage vs \omega')
grid on
legend('Data','Best Fit')
%c

Kt=((Ke/1000)/(pi/30))*141.6; %oz-in/A
ktactual=6.6;%oz-in/A
%% Part 3
File37 = '3.7.xlsx';
data37 = xlsread(File37); %import data
CH037 = data37(5:100004,1); %Ch0
CH137 = data37(5:100004,2); %Ch1
t37 = data37(5:100004,3)-.2935; %time
File39 = '3.9.xlsx';
data39 = xlsread(File39); %import data
CH039 = data39(5:100004,1); %Ch0
CH139 = data39(5:100004,2); %Ch1
t39 = data39(5:100004,3); %time
%a
% system of equations in document

%b

ei=6; %Volts
R=3.6; %Ohms
TStall=(Kt*ei)/R; %ozf-in

%c
Ksys=abs(mean(CH037(end-100:end))/mean(CH137(end-100:end))); %K
Kmot=Ksys/Ktach; %KRPM/Volt for c
tau_Voltage=(.632*(8))-4; 
for i=1:length(t37)
    if CH037(i)  >= tau_Voltage
        tau632_motor = t37(i); %Tau for part c in seconds
        break
       
    else 
    end
end
figure(2)
plot(t37,CH037,t37,CH137,tau632_motor,tau_Voltage,'d')
axis([-.05 .1 -8 8])
grid on
xlabel('Time, t (seconds)')
ylabel('Voltage (Volts)')
legend('Tachometer Output, e_o','Motor Input, e_a')
title('Tachometer Output and Motor Input With Step Change')

%d
R42=4.2;%Ohms
B=.5*(Kt-(Kmot*Ke*Kt))/(Kmot*R42); %oz-in/KRPM
J=.5*((tau632_motor*((Ke*Kt)+(R42*B)))/R42)*(.001/(pi/30));

%e

t39new=t39(79200:end)-.19;
CH039new=CH039(79200:end);
CH139new=CH139(79200:end);
figure(3)
CH039S=wsmooth(CH039new,t39new,6);%Smoothed
CH139S=wsmooth(CH139new,t39new,6);%Smoothed
meanCH039=mean(CH039S(end-100:end));
AMP=(CH039S(1)-CH039S(end));
plot(t39,CH039,t39,CH139)%-meanCH039)%tau632_motor,tau_Voltage,'d')
%plot(t39,CH039)
%axis([-.05 .1 -8 8])
grid on
xlabel('Time, t (seconds)')
ylabel('Voltage (Volts)')
title('Motor Subjected to 5\Omega Resistor')

% for i = 1:length(t39new)
%     
%     if CH039S(i)<= (.368*AMP)+meanCH039
%         tau39 = t39new(i);
%         break
%        
%     else 
%     end
% end
tau39=find((abs((.368*AMP)+meanCH039)-CH039S) <= .001);
plot(t39new,CH039S,t39new(tau39(end)),CH039S(tau39(end)),'d')
grid on
xlabel('Time, t (seconds)')
ylabel('Voltage (Volts)')
title('Motor Subjected to 5\Omega Resistor')
tau39=t39new(tau39(end));
Settling_Time=4*tau39;
sserror=(CH039new(1)-meanCH039)/CH039new(1); %(Initial Voltage - Final Voltage)/(Initial Voltage)
