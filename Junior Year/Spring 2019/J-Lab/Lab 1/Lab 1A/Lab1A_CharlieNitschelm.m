% Lab 1A Code
% There are 4 Main Sections
close all
clear all
%% Part 1

load('Force.mat')
load('Displacement.mat')
Force = unnamed1;
Disp = unnamed;

A = polyfit(Disp,Force,1);


Best_Fit_Force = A(1)*Disp + A(2);

figure(1)
plot(Disp,Force,'o')
hold on
plot(Disp,Best_Fit_Force,'r')
text(1.5,20,strcat('k=' , num2str(A(1),4),' N/mm'),'FontSize',14)
ylabel('Force (N)','FontSize',15)
set(gca,'fontsize',14)
xlabel('Spring Displacement (mm)','FontSize',15)
set(gca,'fontsize',14)
title('Varying Force with Spring Displacement','FontSize',15)
hold off


%% Part 2

Stress_Strain = xlsread('power2.xlsx');

Stress = Stress_Strain(:,2);
Strain = Stress_Strain(:,1);
f = fit(Strain,Stress,'power1');
a = 1527;
b = 0.2368;

figure(2)
plot(f,'--')
hold on
plot(Strain,Stress)
xlabel('\epsilon [ ]','FontSize',16)
set(gca,'fontsize',14)
ylabel('\sigma [MPa]','FontSize',16)
set(gca,'fontsize',14)
text(.15,600,'\sigma_{true} = K\epsilon_{pl}^n')
text(.15,500,strcat('K=',num2str(a,4)))
text(.15,400,strcat('n=',num2str(b,4)))
hold off


%% Part 3

nHeaderLines = 30;
Long = importdata('Long17psi.lvm','\t',nHeaderLines);
time = Long.data(:,1);
voltage = Long.data(:,2);
plot(time,voltage)

for i=1:length(time)
    if time(i) > -0.0035
        basetime = i;
        break
    end
end

baseline = mean(voltage(1:basetime));
basedev  = std(voltage(1:basetime));
threshold = 5*basedev;

for i=1:length(time)
    if (abs(voltage(i) - baseline)>threshold)
        starttime=i;
        break
    end
end

newtime = time(starttime:length(time)) - time(starttime);
newvolt = voltage(starttime:length(time));

smo=20;
mask = ones(smo,1)/smo;
smoothvolt = conv(newvolt,mask,'same');
smoothvolt = smoothvolt(1:end-smo/2);
smoothtime = newtime(1:end-smo/2);

[peakLoc,peakMag] = peakfinder(smoothvolt,.015);


figure (3)
subplot(3,1,1)
plot(newtime,newvolt,'g')
text(0,-0.03,'\leftarrow Start')
xlabel('Time [s]','FontSize',16)
set(gca,'fontsize',14)
ylabel('Voltage [V]','FontSize',16)
set(gca,'fontsize',14)

subplot(3,1,2)
plot(smoothtime,smoothvolt,'--')
text(0,-0.03,'\leftarrow Start')
xlabel('Time [s]','FontSize',16)
set(gca,'fontsize',14)
ylabel('Voltage [V]','FontSize',16)
set(gca,'fontsize',14)

subplot(3,1,3)
plot(smoothtime,smoothvolt,'--')
hold on
plot(smoothtime(peakLoc),peakMag,'x','MarkerSize',14,'MarkerEdgeColor','r')
plot(newtime,newvolt,'g')
text(0,-0.03,'\leftarrow Start')
xlabel('Time [s]','FontSize',16)
set(gca,'fontsize',14)
ylabel('Voltage [V]','FontSize',16)
set(gca,'fontsize',14)
hold off




figure (4)
subplot(2,3,1)
plot(newtime,newvolt,'g',smoothtime,smoothvolt,'--')
hold on
plot(smoothtime(peakLoc),peakMag,'x','MarkerSize',14,'MarkerEdgeColor','r')
hold off
ylim([1.52,1.65])
xlim([0.015,0.025])
set(gca,'fontsize',14)
ylabel('Voltage [V]','FontSize',16)
set(gca,'fontsize',14)

subplot(2,3,2)
plot(newtime,newvolt,'g',smoothtime,smoothvolt,'--')
hold on
plot(smoothtime(peakLoc),peakMag,'x','MarkerSize',14,'MarkerEdgeColor','r')
hold off
ylim([1.52,1.65])
xlim([.028,.038])
set(gca,'fontsize',14)
set(gca,'fontsize',14)

subplot(2,3,3)
plot(newtime,newvolt,'g',smoothtime,smoothvolt,'--')
hold on
plot(smoothtime(peakLoc),peakMag,'x','MarkerSize',14,'MarkerEdgeColor','r')
hold off
ylim([1.52,1.65])
xlim([.042,.052])
set(gca,'fontsize',14)
set(gca,'fontsize',14)

subplot(2,3,4)
plot(newtime,newvolt,'g',smoothtime,smoothvolt,'--')
hold on
plot(smoothtime(peakLoc),peakMag,'x','MarkerSize',14,'MarkerEdgeColor','r')
hold off
ylim([1.52,1.65])
xlim([.055,.065])
xlabel('Time [s]','FontSize',16)
set(gca,'fontsize',14)
ylabel('Voltage [V]','FontSize',16)
set(gca,'fontsize',14)

subplot(2,3,5)
plot(newtime,newvolt,'g',smoothtime,smoothvolt,'--')
hold on
plot(smoothtime(peakLoc),peakMag,'x','MarkerSize',14,'MarkerEdgeColor','r')
hold off
ylim([1.52,1.65])
xlim([.07,.08])
xlabel('Time [s]','FontSize',16)
set(gca,'fontsize',14)
set(gca,'fontsize',14)

subplot(2,3,6)
plot(newtime,newvolt,'g',smoothtime,smoothvolt,'--')
hold on
plot(smoothtime(peakLoc),peakMag,'x','MarkerSize',14,'MarkerEdgeColor','r')
hold off
ylim([1.52,1.65])
xlim([.083,.093])
xlabel('Time [s]','FontSize',16)
set(gca,'fontsize',14)
set(gca,'fontsize',14)


tablehelper = newtime(peakLoc);
f = figure('Position',[200,200,400,150]);
dat = [tablehelper(1) peakMag(1); tablehelper(2) peakMag(2); tablehelper(3) peakMag(3); tablehelper(4) peakMag(4); tablehelper(5) peakMag(5); tablehelper(6) peakMag(6)]; 
cnames = {'Time (s)', 'Amplitude (V)'};
rnames = {'First Peak','Second Peak','Third Peak','Fourth Peak','Fifth Peak','Sixth Peak'};
t = uitable('Parent',f,'Data',dat,'ColumnName',cnames,'RowName',rnames,'Position',[20 20 375 135]);




