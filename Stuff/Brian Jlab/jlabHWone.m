close all; clear all; clc;

x = xlsread('jlabHWone.xlsx','A2:A19');
y = xlsread('jlabHWone.xlsx','B2:B19');

StandDev = std(y);
p = polyfit(x,y,1);
p1 = p(1);
p2 = p(2);
yFIT = p(1)*x + p(2);

t = 2.11;
sampleNumber = 18;
Sx = StandDev / sqrt(sampleNumber);
conf = t*Sx;
StandDevFit  = std(yFIT);
SxFIT = StandDevFit / sqrt(sampleNumber);
confFIT = t*SxFIT;

DataConfP = y + conf;
DataConfM = y - conf;
BFConfP   = yFIT + confFIT;
BFConfM   = yFIT - confFIT;

plot(x,y,'o'); hold on;
plot(x,yFIT,'--'); hold on;
% plot(x,DataConfP); hold on;
% plot(x,DataConfM); hold on;
% plot(x,BFConfP); hold on;
% plot(x,BFConfM); hold on;
grid minor
% legend('Data Points','Best Fit Line','Upper Data Confidence Limit','Location','Northwest')

xGuess = 45;
yFITGuess = (p(1)*xGuess + p(2));
GuessHigh = yFITGuess + confFIT;
GusesLow  = yFITGuess - confFIT;




P = polyfit(x,y,4);
Sum = 0;
for i = 1:length(x)
    Sum(i) =(x(i)-((P(1)*x(i)^4) + (P(2)*x(i)^3) + (P(3)*x(i)^2) + (P(4)*x(i)) + P(5)))^2;
end
Total = sum(Sum);

Order = 4;
v = sampleNumber - (Order+1);
Sxy = sqrt(Total/v);

sum1 = 0;
for i = 1:length(x)
    sum1(i) = (x(i) - mean(x))^2;
end
Total1 = sum(sum1);

for i = 1:length(x)
    Yc(i) = P(5) + (P(4)*x(i)) + (P(3)*x(i)^2) + (P(2)*x(i)^3) + (P(1)*x(i)^4);
end

for i = 1:length(x)
YxP(i) = Yc(i) + ((t*Sxy)*((1/sampleNumber)+((x(i) - mean(x))^2/Total1))^(1/2));
YxM(i) = Yc(i) - ((t*Sxy)*((1/sampleNumber)+((x(i) - mean(x))^2/Total1))^(1/2));
end
plot(x,YxP); hold on;
plot(x,YxM); hold on;


for i = 1:length(x)
YxPF(i) = Yc(i) + ((t*Sxy)*(1+((1/sampleNumber)+((x(i) - mean(x))^2)/Total1))^(1/2));
YxMF(i) = Yc(i) - ((t*Sxy)*(1+((1/sampleNumber)+((x(i) - mean(x))^2)/Total1))^(1/2));
end

plot(x,YxPF); hold on;
plot(x,YxMF);






