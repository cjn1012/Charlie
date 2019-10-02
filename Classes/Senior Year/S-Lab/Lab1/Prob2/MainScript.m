clear all
close all
d = abs(-.140723854 - 0.140276146);
e = abs(-0.140724592 - 0.140275408);
Time_2 = linspace(0,d,100000);
Volt_in2 = xlsread('2Volt.csv','2Volt' ,'A8:A100007'); % 2Volt
Volt_out2 = xlsread('2Volt.csv','2Volt' ,'B8:B100007'); % 2Volt out
Time_3 = linspace(0,e,100000);
Volt_in3 = xlsread('-3-5Volt.csv','-3-5Volt' ,'A8:A100007'); % 3Volt
Volt_out3 = xlsread('-3-5Volt.csv','-3-5Volt' ,'B8:B100007'); % 3Volt out

%% +- 2 volts
Volt2 = Volt_in2(14501:20000);
Volt2out = Volt_out2(14501:20000);
Time2 = Time_2(1:5500);

% Inserting peaks and calc constants
Constant=0.05;
[Locations,Values]=peakfinder(Volt2out,Constant);
Value_Log = Values-2;
Log = log(Value_Log(1)./Value_Log(1:2));
Con = transpose(1:2);
[mb,Time_Wave] = polyfit(Con,Log,1);
Best_Fit_3 = mb(1)*Con+mb(2);
alpha = mb(1);
Damp_Ratio2 = (alpha/sqrt(4*pi^2+alpha^2));

Period = (Time2(Locations(3))-Time2(Locations(1)))/2;
Damping_Ratio_2 = 2*pi/Period;
Damp_Natural_Freq_2 = mean(Damping_Ratio_2);
Undamped_Natural_Freq_2 = Damp_Natural_Freq_2/(sqrt(1-Damp_Ratio2^2));


% figure(1)
% plot(Time2,Volt2,Time2,Volt2out,Time2(Locations),Values,'o')


% Part c

R2 = 20150; %ohms
L2 =((2*Damp_Ratio2*R2)/Undamped_Natural_Freq_2); %H
C2 = 1/(L2*((Undamped_Natural_Freq_2)^2)); % Micro Coulumb

% Part d

num2 = [R2];
den2 = [L2*C2*R2 L2 R2];
sim('Simu20')
Theo_Time2 = ScopeData(:,1);
Theo_Volt2 = ScopeData(:,2)-2;

figure(2)
plot(Time2,Volt2out,Theo_Time2,Theo_Volt2,Time2(Locations),Values,'o')
title('Step Response of a Second Order System for -2 to +2 Volts')
xlabel('Time (s)')
ylabel('Volts (V)')
legend('Experimental Data','Theoretical Data')

%% -3 to 5 volt

Volt3 = Volt_in3(14501:20000);
Volt3out = Volt_out3(14501:20000);
Time3 = Time_3(1:5500);

% Inserting peaks and calc constants
Constant=0.05;
[Locations,Values]=peakfinder(Volt3out,Constant);
Value_Log = Values-5;
Log = log(Value_Log(1)./Value_Log(1:2));
Con = transpose(1:2);
[mb,Time_Wave] = polyfit(Con,Log,1);
Best_Fit_3 = mb(1)*Con+mb(2);
alpha = mb(1);
Damp_Ratio3 = (alpha/sqrt(4*pi^2+alpha^2));

Period = (Time3(Locations(4))-Time3(Locations(1)))/3;
Damping_Ratio_3 = 2*pi/Period;
Damp_Natural_Freq_3 = mean(Damping_Ratio_3);
Undamped_Natural_Freq_3 = Damp_Natural_Freq_3/(sqrt(1-Damp_Ratio3^2));


% figure(3)
% plot(Time3,Volt3,Time3,Volt3out,Time3(Locations),Values,'o')


% Part c

R3 = 20150; %ohms
L3 =((2*Damp_Ratio3*R3)/Undamped_Natural_Freq_3); %H
C3 = 1/(L3*((Undamped_Natural_Freq_3)^2)); % Micro Coulumb

% Part d

num3 = [R3];
den3 = [L3*C3*R3 L3 R3];
sim('Simu3')
Theo_Time3 = ScopeData(:,1);
Theo_Volt3 = ScopeData(:,2)-3;

figure(4)
plot(Time3,Volt3out,Theo_Time3,Theo_Volt3,Time3(Locations),Values,'o')
title('Step Response of a Second Order System for -3 to 5 Volts')
xlabel('Time (s)')
ylabel('Volts (V)')
legend('Experimental Data','Theoretical Data')




%% Part 2 of 2 Freqiency response

% +-2v

% Experimental Data
Mag = [0.021796472,0.341504021,1.981738645,9.502997475,-2.021055162,-8.453984954,-21.04216733];
Frq = [14.25,71.25,142.5,285,427.5,570,1140]*2*pi;
PhA = [-0.25651988,-1.542857143,-6.942857143,-72,-156.5217391,-172.8,-172.4059293];


sys2 = tf(num2,den2);

[Magnitude,Phase,Frequency] = bode(sys2);
Phase = Phase(:);
Magnitude = Magnitude(:);
Magnitude = 20*log10(abs(Magnitude));
figure(5)
subplot(2,1,1)
semilogx(Frequency,Magnitude(:),Frq,Mag,'-o')
xlim([70,10000])
ylabel('Magnitude (dB)')

subplot(2,1,2)
semilogx(Frequency,Phase,Frq,PhA,'-o')
xlim([70,10000])
legend('Theoretical Data', 'Experimental Data')
xlabel('Frequency (rad/s)')
ylabel('Phase Angle (deg)')


% part c

for x = 1:length(Phase)
    if Phase(x)<-90
        nat_freq = (Frequency(x-1)+ ((-90-Phase(x-1))/(Phase(x)-Phase(x-1)))*(Frequency(x)-Frequency(x-1)));
        break
    end 
end

for x = 2:length(Mag(1:6))
    r = Frq(x)/nat_freq;
    damp(x-1) = (tand(-PhA(x))*(1-r^2))/2*r;
end
damping = mean(damp)

R2_1 = 20150; %ohms
L2_1 =((2*damping*R2_1)/nat_freq); %H
C2_1 = 1/(L2_1*((nat_freq)^2)); % Micro Coulumb

num4 = [R2_1];
den4 = [L2_1*C2_1*R2_1 L2_1 R2_1];
sys4 = tf(num4,den4);

[Magnitude,Phase,Frequency] = bode(sys4);
Phase = Phase(:);
Magnitude = Magnitude(:);
Magnitude = 20*log10(abs(Magnitude));

%% Part 3 second order

Frq2 = xlsread('Frequency Response Data.xlsx','2nd Order Magnitude' ,'A12:A1013');
Frq2 = Frq2.*2*pi
Mag2 = xlsread('Frequency Response Data.xlsx','2nd Order Magnitude' ,'B12:B1013'); % 2Volt
Pha2 = xlsread('Frequency Response Data.xlsx','2nd Order Phase' ,'B12:B1013'); % 2Volt

figure(6)
subplot(2,1,1)
semilogx(Frequency,Magnitude(:),Frq2,Mag2)
xlim([70,10000])
ylabel('Magnitude (dB)')

subplot(2,1,2)
semilogx(Frequency,Phase,Frq2,Pha2)
xlim([70,10000])
legend('Theoretical Data', 'Experimental Data')
xlabel('Frequency (rad/s)')
ylabel('Phase Angle (deg)')

for x = 1:length(Pha2)
    if Pha2(x)<-90
        nat_freq = (Frq2(x-1)+ ((-90-Pha2(x-1))/(Pha2(x)-Pha2(x-1)))*(Frq2(x)-Frq2(x-1)));
        break
    end 
end

for x = 2:length(Mag2(1:150))
    r = Frq2(x)/nat_freq;
    damp2(x-1) = (tand(-Pha2(x))*(1-r^2))/2*r;
end
damping2 = mean(damp2)


R2_2 = 20150; %ohms
L2_2 =((2*damping2*R2_2)/nat_freq); %H
C2_2 = 1/(L2_2*((nat_freq)^2)); % Micro Coulumb
