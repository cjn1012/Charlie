clear all
close all
clc

%% Question 1
%plot a 1st order response with a time constant of 0.5sec
tau1=0.5;

t=linspace(0,5,100);

T_final=273.15;
T_start=76+273.15;
% Temp_tau_std_value_BBI(i)=BBI_temp_final-((BBI_temp_final-Temp_start_poly_BBI)*exp(-new_time_start_std_BBI(i)/Tau_std_value_BBI));
for i=1:length(t)
    temp1(i)=T_final-((T_final-T_start)*exp(-t(i)/tau1));
end
T_5deg=5+273.15;
time_5deg=-tau1*log((T_5deg-T_final)/(T_start-T_final));
plot(t,temp1)
text(3, 340, ['t(278.15K)= ' num2str(time_5deg) ' sec'])
% ylim([.5, 80])
grid minor
title('Time vs Temperature Plot (\tau = .5sec)')
xlabel('Time (sec)')
ylabel('Temperature (K)')

%% Question 4
tau4=0.05;
w4=linspace(0,5000,10000);
for i=1:length(w4)
   m4(i)=1/(sqrt(1+(w4(i)*tau4)^2));
end
figure
plot(w4,m4);
title('\omega vs M(\omega)')
xlabel('\omega (rad/sec)')
ylabel('M(\omega)')
grid minor

 