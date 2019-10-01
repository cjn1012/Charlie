clear all
close all

Time_2 = linspace(-.140723854,0.140276146,100000);
Volt_in2 = xlsread('2Volt.csv','2Volt' ,'A8:A100007'); % 2Volt
Volt_out2 = xlsread('2Volt.csv','2Volt' ,'B8:B100007'); % 2Volt
Time_3 = linspace(-0.140724592,0.140275408,100000);
Volt_in3 = xlsread('-3-5Volt.csv','-3-5Volt' ,'A8:A100007'); % 3Volt
Volt_out3 = xlsread('-3-5Volt.csv','-3-5Volt' ,'B8:B100007'); % 3Volt

%% +- 2 volts
Volt2 = Volt_in2(14494:20000);
Volt2out = Volt_out2(14494:20000);
Time2 = Time_2(14494:20000);

Constant=0.03;
[Locations,Values]=peakfinder(Volt2,Constant);
Locations = Locations(2:end);
Values = Values(2:end);

Log = log(Values(2)./Values(2:length(Values)));
Con = transpose(2:length(Values));
[mb,Time_Wave] = polyfit(Con,Log,1);
Best_Fit_3 = mb(1)*Con+mb(2);
alpha = mb(1);
Damp_Ratio2 = (alpha/sqrt(4*pi^2+alpha^2));

Damping_Ratio_2 = zeros(length(Values),1);
for n = 1:length(Values)-1
    Period1 = Time2(Locations(n+1))-Time2(Locations(n));
    Damping_Ratio_2(n) = (2*pi)/Period1;
end
Damp_Natural_Freq_2 = mean(Damping_Ratio_2);
Undamped_Natural_Freq_2 = Damp_Natural_Freq_2/(sqrt(1-Damp_Ratio2^2));



figure(1)
plot(Time2,Volt2,Time2,Volt2out,Time2(Locations),Values,'o')
figure(2)
plot(Time_2,Volt_in2,Time_2,Volt_out2)


%% -3 +5 volts

Volt3 = Volt_in3(14501:20000);
Time3 = Time_3(14501:20000);

Constant=0.03;
[Locations,Values]=peakfinder(Volt3,Constant);
Locations = Locations(2:end);
Values = Values(2:end);

Log = log(Values(2)./Values(2:length(Values)));
Con = transpose(2:length(Values));
[mb,Time_Wave] = polyfit(Con,Log,1);
Best_Fit_3 = mb(1)*Con+mb(2);
alpha = mb(1);
Damp_Ratio3 = (alpha/sqrt(4*pi^2+alpha^2));

Damping_Ratio_3 = zeros(length(Values),1);
for n = 1:length(Values)-1
    Period3 = Time3(Locations(n+1))-Time3(Locations(n));
    Damping_Ratio_3(n) = (2*pi)/Period3;
end
Damp_Natural_Freq_3 = mean(Damping_Ratio_3);
Undamped_Natural_Freq_3 = Damp_Natural_Freq_3/(sqrt(1-Damp_Ratio3^2));

figure(2)
plot(Time3,Volt3,Time3(Locations),Values,'o')
















%% Part c

R2 = 20150 %ohms
L2 =(2*Damp_Ratio2*R2)/Undamped_Natural_Freq_2; %H
C2 = 1/(L2*(Undamped_Natural_Freq_2^2))%*10^6 % Micro Coulumb

R3 = 20150 %ohms
L3 =(2*Damp_Ratio3*R3)/Undamped_Natural_Freq_3; %H
C3 = 1/(L3*(Undamped_Natural_Freq_3^2))%*10^6 % Micro Coulumb


%% Part d

num2 = [R2*C2 1]
den2 = [L2*C2 R2*C2 1]
sim('Simu2')
Theo_Time2 = ScopeData(:,1)
Theo_Volt2 = ScopeData(:,2)


num3 = [R3*C3 1]
den3 = [L3*C3 R3*C3 1]
sim('Simu3')
Theo_Time3 = ScopeData(:,1)
Theo_Volt3 = ScopeData(:,2)


figure(3)
plot(Time2,Volt2,Theo_Time2,Theo_Volt2)



















