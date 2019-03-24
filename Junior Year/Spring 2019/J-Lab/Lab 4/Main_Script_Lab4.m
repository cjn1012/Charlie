% Lab 4 %
clear all
close all

%% Import Data

% Reading in entire data set
Header=29;
Small_Flood    = importdata('SmallPipe_Flood.lvm','\t',Header);
Small_Balloon  = importdata('SmallPipe_Balloon.lvm','\t',Header);
Medium_Flood   = importdata('MediumPipe_Flood.lvm','\t',Header);
Medium_Balloon = importdata('MediumPipe_Balloon.lvm','\t',Header);
Long_Flood     = importdata('LongPipe_Flood.lvm','\t',Header);
Long_Balloon   = importdata('LongPipe_Balloon.lvm','\t',Header);

% Separating entire data sets into vectors of Time and Voltage
Time_Small_Flood    = Small_Flood.data(:,1);
Time_Small_Balloon  = Small_Balloon.data(:,1);
Time_Medium_Flood   = Medium_Flood.data(:,1);
Time_Medium_Balloon = Medium_Balloon.data(:,1);
Time_Long_Flood     = Long_Flood.data(:,1);
Time_Long_Balloon   = Long_Balloon.data(:,1);

Voltage_Small_Flood    = Small_Flood.data(:,2);
Voltage_Small_Balloon  = Small_Balloon.data(:,2);
Voltage_Medium_Flood   = Medium_Flood.data(:,2);
Voltage_Medium_Balloon = Medium_Balloon.data(:,2);
Voltage_Long_Flood     = Long_Flood.data(:,2);
Voltage_Long_Balloon   = Long_Balloon.data(:,2);

%% Flood Tests

% Small Pipe

T = Time_Small_Flood(2)-Time_Small_Flood(1); % Calculated time interval between data points
Freq = 1/T; % Sampling Frequency
Length = size(Voltage_Small_Flood); % Number of points in vector
Power_2 = 2^nextpow2(Length(1)); % power of 2 from length y
FFT = fft(Voltage_Small_Flood,Power_2)./Length(1); % fft bull
Spaced_Points = Freq/2*linspace(0,1,Power_2/2+1); % spaced point vector using length
FFT_A = 20*log10(abs(FFT(1:Power_2/2+1))); % amplitude
Length_2 = length(Voltage_Small_Flood);
Freq_Reso_Flood_Small = 1/(Length_2*T); %FOR GRAPH

a = 343;
l = .17;
V = l*pi*.0042^2;
Vt = 6.5*10^-8;
Damp_Natural_Freq_Small_Flood = a/(l*sqrt(0.5+(V/Vt))); % FOR GRAPH

figure(1)

subplot(2,1,1)
plot (Time_Small_Flood,Voltage_Small_Flood);
title('Short Tube Flood Valve Response')
xlabel ('Time (sec)')
ylabel ('Volts (V)')
text(.03,2,strcat('Resolution =',num2str(Freq_Reso_Flood_Small),'1/s'))
text(.03,1.5,strcat('Damped Natural Frequency = ',num2str(Damp_Natural_Freq_Small_Flood),'Hz'))
xlim([0,.1])

subplot(2,1,2)
semilogx(Spaced_Points,FFT_A);
title('Single Sided Amplitude Spectrum')
xlabel ('Frequency (Hz)')
ylabel ('Log Magnitude (dBv)')


% Medium Pipe

T = Time_Medium_Flood(2)-Time_Medium_Flood(1); % Calculated time interval between data points
Freq = 1/T; % Sampling Frequency
Length = size(Voltage_Medium_Flood); % Number of points in vector
Power_2 = 2^nextpow2(Length(1)); % power of 2 from length y
FFT = fft(Voltage_Medium_Flood,Power_2)./Length(1); % fft bull
Spaced_Points = Freq/2*linspace(0,1,Power_2/2+1); % spaced point vector using length
FFT_A = 20*log10(abs(FFT(1:Power_2/2+1))); % amplitude
Length_2 = length(Voltage_Medium_Flood);
Freq_Reso_Flood_Medium = 1/(Length_2*T); %FOR GRAPH

a = 343;
l = .795;
V = l*pi*.0042^2;
Vt = 6.5*10^-8;
Damp_Natural_Freq_Medium_Flood = a/(l*sqrt(0.5+(V/Vt)));

figure(2)

subplot(2,1,1)
plot (Time_Medium_Flood,Voltage_Medium_Flood);
title('Medium Tube Flood Valve Response')
xlabel ('Time (sec)')
ylabel ('Volts (V)')
text(.03,2,strcat('Resolution =',num2str(Freq_Reso_Flood_Medium),'1/s'))
text(.03,1.5,strcat('Damped Natural Frequency = ',num2str(Damp_Natural_Freq_Medium_Flood),'Hz'))
xlim([0,.1])

subplot(2,1,2)
semilogx(Spaced_Points,FFT_A);
title('Single Sided Amplitude Spectrum')
xlabel ('Frequency (Hz)')
ylabel ('Log Magnitude (dBv)')



% Long Pipe

T = Time_Long_Flood(2)-Time_Long_Flood(1); % Calculated time interval between data points
Freq = 1/T; % Sampling Frequency
Length = size(Voltage_Small_Flood); % Number of points in vector
Power_2 = 2^nextpow2(Length(1)); % power of 2 from length y
FFT = fft(Voltage_Small_Flood,Power_2)./Length(1); % fft bull
Spaced_Points = Freq/2*linspace(0,1,Power_2/2+1); % spaced point vector using length
FFT_A = 20*log10(abs(FFT(1:Power_2/2+1))); % amplitude
Length_2 = length(Voltage_Long_Flood);
Freq_Reso_Flood_Long = 1/(Length_2*T); %FOR GRAPH

a = 343;
l = 1.05;
V = l*pi*.0042^2;
Vt = 6.5*10^-8;
Damp_Natural_Freq_Long_Flood = a/(l*sqrt(0.5+(V/Vt)));

figure(3)

subplot(2,1,1)
plot (Time_Long_Flood,Voltage_Long_Flood);
title('Long Tube Flood Valve Response')
xlabel ('Time (sec)')
ylabel ('Volts (V)')
text(.03,2,strcat('Resolution =',num2str(Freq_Reso_Flood_Long),'1/s'))
text(.03,1.5,strcat('Damped Natural Frequency = ',num2str(Damp_Natural_Freq_Long_Flood),'Hz'))
xlim([0,.1])

subplot(2,1,2)
semilogx(Spaced_Points,FFT_A);
title('Single Sided Amplitude Spectrum')
xlabel ('Frequency (Hz)')
ylabel ('Log Magnitude (dBv)')


%% Balloon Tests

% Small Pipe

T = Time_Small_Balloon(2)-Time_Small_Balloon(1); % Calculated time interval between data points
Freq = 1/T; % Sampling Frequency
Length = size(Voltage_Small_Balloon); % Number of points in vector
Power_2 = 2^nextpow2(Length(1)); % power of 2 from length y
FFT = fft(Voltage_Small_Balloon,Power_2)./Length(1); % fft bull
Spaced_Points = Freq/2*linspace(0,1,Power_2/2+1); % spaced point vector using length
FFT_A = 20*log10(abs(FFT(1:Power_2/2+1))); % amplitude
Length_2 = length(Voltage_Small_Balloon);
Freq_Reso_Balloon_Small = 1/(Length_2*T); %FOR GRAPH

a = 343;
l = .17;
V = l*pi*.0042^2;
Vt = 6.5*10^-8;
Damp_Natural_Freq_Small_Balloon = a/(l*sqrt(0.5+(V/Vt)));

figure(4)

subplot(2,1,1)
plot (Time_Small_Balloon,Voltage_Small_Balloon);
title('Short Tube Balloon Response')
xlabel ('Time (sec)')
ylabel ('Volts (V)')
text(.015,-.03,strcat('Resolution =',num2str(Freq_Reso_Balloon_Small),'1/s'))
text(.015,-.06,strcat('Damped Natural Frequency = ',num2str(Damp_Natural_Freq_Small_Balloon),'Hz'))
xlim([-.025,.075])

subplot(2,1,2)
semilogx(Spaced_Points,FFT_A);
title('Single Sided Amplitude Spectrum')
xlabel ('Frequency (Hz)')
ylabel ('Log Magnitude (dBv)')


% Medium Pipe

T = Time_Medium_Balloon(2)-Time_Medium_Balloon(1); % Calculated time interval between data points
Freq = 1/T; % Sampling Frequency
Length = size(Voltage_Medium_Balloon); % Number of points in vector
Power_2 = 2^nextpow2(Length(1)); % power of 2 from length y
FFT = fft(Voltage_Medium_Balloon,Power_2)./Length(1); % fft bull
Spaced_Points = Freq/2*linspace(0,1,Power_2/2+1); % spaced point vector using length
FFT_A = 20*log10(abs(FFT(1:Power_2/2+1))); % amplitude
Length_2 = length(Voltage_Medium_Balloon);
Freq_Reso_Balloon_Medium = 1/(Length_2*T); %FOR GRAPH

a = 343;
l = .795;
V = l*pi*.0042^2;
Vt = 6.5*10^-8;
Damp_Natural_Freq_Medium_Balloon = a/(l*sqrt(0.5+(V/Vt)));

figure(5)

subplot(2,1,1)
plot (Time_Medium_Balloon,Voltage_Medium_Balloon);
title('Medium Tube Balloon Response')
xlabel ('Time (sec)')
ylabel ('Volts (V)')
text(.015,-.07,strcat('Resolution =',num2str(Freq_Reso_Balloon_Medium),'1/s'))
text(.015,-.10,strcat('Damped Natural Frequency = ',num2str(Damp_Natural_Freq_Medium_Balloon),'Hz'))
xlim([-.025,.075])

subplot(2,1,2)
semilogx(Spaced_Points,FFT_A);
title('Single Sided Amplitude Spectrum')
xlabel ('Frequency (Hz)')
ylabel ('Log Magnitude (dBv)')



% Long Pipe

T = Time_Long_Balloon(2)-Time_Long_Balloon(1); % Calculated time interval between data points
Freq = 1/T; % Sampling Frequency
Length = size(Voltage_Small_Balloon); % Number of points in vector
Power_2 = 2^nextpow2(Length(1)); % power of 2 from length y
FFT = fft(Voltage_Small_Balloon,Power_2)./Length(1); % fft bull
Spaced_Points = Freq/2*linspace(0,1,Power_2/2+1); % spaced point vector using length
FFT_A = 20*log10(abs(FFT(1:Power_2/2+1))); % amplitude
Length_2 = length(Voltage_Long_Balloon);
Freq_Reso_Balloon_Long = 1/(Length_2*T); %FOR GRAPH

a = 343;
l = 1.05;
V = l*pi*.0042^2;
Vt = 6.5*10^-8;
Damp_Natural_Freq_Long_Balloon = a/(l*sqrt(0.5+(V/Vt)));

figure(6)

subplot(2,1,1)
plot (Time_Long_Balloon,Voltage_Small_Balloon);
title('Long Tube Balloon Response')
xlabel ('Time (sec)')
ylabel ('Volts (V)')
text(.015,-.05,strcat('Resolution =',num2str(Freq_Reso_Balloon_Long),'1/s'))
text(.015,-.09,strcat('Damped Natural Frequency = ',num2str(Damp_Natural_Freq_Long_Balloon),'Hz'))
xlim([-.025,.075])

subplot(2,1,2)
semilogx(Spaced_Points,FFT_A);
title('Single Sided Amplitude Spectrum')
xlabel ('Frequency (Hz)')
ylabel ('Log Magnitude (dBv)')








%% Peakfinder

% Small Pipe
Constant = 0.028;
[Volt_Locations_Small,Volt_Small]=peakfinder(Voltage_Small_Flood,Constant);
Damping_Ratio_Small_Flood = zeros(length(Volt_Small),1);
for n = 1:length(Volt_Small)-1
    Period_Small = Time_Small_Flood(Volt_Locations_Small(n+1))-Time_Small_Flood(Volt_Locations_Small(n));
    Damping_Ratio_Small_Flood(n) = 2*pi/Period_Small;
end
Damp_Mean_Small = mean(Damping_Ratio_Small_Flood);
figure(7)
plot(Time_Small_Flood,Voltage_Small_Flood,Time_Small_Flood(Volt_Locations_Small),Volt_Small,'d')
title('Short Tube Flood Valve Response with PeakFinder')
text(.04,2.35,strcat('Damped Natural Frequency =',num2str(Damp_Mean_Small),'Hz'))
xlabel ('Time (sec)')
ylabel ('Volts (V)')
xlim([.008,.15])


% Medium Pipe
Constant = 0.036;
[Volt_Locations_Medium,Volt_Medium]=peakfinder(Voltage_Medium_Flood,Constant);
Damping_Ratio_Medium_Flood = zeros(length(Volt_Medium),1);
for n = 1:length(Volt_Medium)-1
    Period_Medium = Time_Medium_Flood(Volt_Locations_Medium(n+1))-Time_Medium_Flood(Volt_Locations_Medium(n));
    Damping_Ratio_Medium_Flood(n) = 2*pi/Period_Medium;
end
Damp_Mean_Medium = mean(Damping_Ratio_Medium_Flood);
figure(8)
plot(Time_Medium_Flood,Voltage_Medium_Flood,Time_Medium_Flood(Volt_Locations_Medium),Volt_Medium,'d')
title('Medium Tube Flood Valve Response with PeakFinder')
text(.04,2.35,strcat('Damped Natural Frequency =',num2str(Damp_Mean_Medium),'Hz'))
xlabel ('Time (sec)')
ylabel ('Volts (V)')
xlim([.01,.15])

% Long Pipe
Constant = 0.035;
[Volt_Locations_Long,Volt_Long]=peakfinder(Voltage_Long_Flood,Constant);
Damping_Ratio_Long_Flood = zeros(length(Volt_Long),1);
for n = 1:length(Volt_Long)-1
    Period_Long = Time_Long_Flood(Volt_Locations_Long(n+1))-Time_Long_Flood(Volt_Locations_Long(n));
    Damping_Ratio_Long_Flood(n) = 2*pi/Period_Long;
end
Damp_Mean_Long = mean(Damping_Ratio_Long_Flood);
figure(9)
plot(Time_Long_Flood,Voltage_Long_Flood,Time_Long_Flood(Volt_Locations_Long),Volt_Long,'d')
title('Long Tube Flood Valve Response with PeakFinder')
text(.04,2.35,strcat('Damped Natural Frequency =',num2str(Damp_Mean_Long),'Hz'))
xlabel ('Time (sec)')
ylabel ('Volts (V)')
xlim([.01,.15])



%% p-t response and predicted

Volt_Locations_Small;
Volt_Locations_Medium;
Volt_Locations_Long;
Volt_Small;
Volt_Medium;
Volt_Long;



Constant = 0.01;
[Volt_Locations_Balloon_Small,Volt_Small_Balloon]=peakfinder(Voltage_Small_Balloon,Constant);
Volt_Locations_Balloon_Small = Volt_Locations_Balloon_Small(2:end-1);
Volt_Small_Balloon = Volt_Small_Balloon(2:end-1);
figure(10)
plot(Time_Small_Balloon,Voltage_Small_Balloon,Time_Small_Balloon(Volt_Locations_Balloon_Small),Volt_Small_Balloon,'d')

Constant = 0.02;
[Volt_Locations_Balloon_Medium,Volt_Medium_Balloon]=peakfinder(Voltage_Medium_Balloon,Constant);
Volt_Locations_Balloon_Medium = Volt_Locations_Balloon_Medium(3:end);
Volt_Medium_Balloon = Volt_Medium_Balloon(3:end);
figure(11)
plot(Time_Medium_Balloon,Voltage_Medium_Balloon,Time_Medium_Balloon(Volt_Locations_Balloon_Medium),Volt_Medium_Balloon,'d')

Constant = 0.02;
[Volt_Locations_Balloon_Long,Volt_Long_Balloon]=peakfinder(Voltage_Long_Balloon,Constant);
Volt_Locations_Balloon_Long = Volt_Locations_Balloon_Long(3:end);
Volt_Long_Balloon = Volt_Long_Balloon(3:end);
figure(12)
plot(Time_Long_Balloon,Voltage_Long_Balloon,Time_Long_Balloon(Volt_Locations_Balloon_Long),Volt_Long_Balloon,'d')

Damp_Small_Valve = zeros(length(Volt_Locations_Small)-1,1);
Damp_Medium_Valve = zeros(length(Volt_Locations_Medium)-1,1);
Damp_Long_Valve = zeros(length(Volt_Locations_Long)-1,1);
Damp_Small_Balloon = zeros(length(Volt_Locations_Balloon_Small)-1,1);
Damp_Medium_Balloon = zeros(length(Volt_Locations_Balloon_Medium)-1,1);
Damp_Long_Balloon = zeros(length(Volt_Locations_Balloon_Long)-1,1);

for n=2:length(Volt_Locations_Small)
    Num = (1/(n-1))*log(Volt_Small(2)./Volt_Small(n));
    Damp_Small_Valve(n-1) = Num./sqrt(4*pi^2+Num^2);
end

for n=2:length(Volt_Locations_Medium)
    Num = (1/(n-1))*log(Volt_Medium(2)./Volt_Medium(n));
    Damp_Medium_Valve(n-1) = Num./sqrt(4*pi^2+Num^2);
end

for n=2:length(Volt_Locations_Long)
    Num = (1/(n-1))*log(Volt_Long(2)./Volt_Long(n));
    Damp_Long_Valve(n-1) = Num./sqrt(4*pi^2+Num^2);
end

for n=2:length(Volt_Locations_Balloon_Small)
    Num = (1/(n-1))*abs(log(Volt_Small_Balloon(2)./Volt_Small_Balloon(n)));
    Damp_Small_Balloon(n-1) = Num./sqrt(4*pi^2+Num^2);
end

for n=2:length(Volt_Locations_Balloon_Medium)
    Num = (1/(n-1))*log(Volt_Medium_Balloon(2)./Volt_Medium_Balloon(n));
    Damp_Medium_Balloon(n-1) = Num./sqrt(4*pi^2+Num^2);
end

for n=2:length(Volt_Locations_Balloon_Long)
    Num = (1/(n-1))*abs(log(Volt_Long_Balloon(2)./Volt_Long_Balloon(n)))
    Damp_Long_Balloon(n-1) = Num./sqrt(4*pi^2+Num^2);
end



Ratio_Valve_Small = mean(Damp_Small_Valve);
Ratio_Valve_Medium = mean(Damp_Medium_Valve);
Ratio_Valve_Long = mean(Damp_Long_Valve);
Ratio_Balloon_Small = mean(Damp_Small_Balloon);
Ratio_Balloon_Medium = mean(Damp_Medium_Balloon);
Ratio_Balloon_Long = mean(Damp_Long_Balloon);


Damp_Mean_Small;
Damp_Mean_Medium;
Damp_Mean_Long;

for n = 1:length(Volt_Small_Balloon)-1
    Period_Small = Time_Small_Balloon(Volt_Locations_Balloon_Small(n+1))-Time_Small_Balloon(Volt_Locations_Balloon_Small(n));
    Damping_Ratio_Small_Balloon(n) = 2*pi/Period_Small;
end
Damp_Mean_Small_Balloon = mean(Damping_Ratio_Small_Balloon);

for n = 1:length(Volt_Medium_Balloon)-1
    Period_Medium = Time_Medium_Balloon(Volt_Locations_Balloon_Medium(n+1))-Time_Medium_Balloon(Volt_Locations_Balloon_Medium(n));
    Damping_Ratio_Medium_Balloon(n) = 2*pi/Period_Medium;
end
Damp_Mean_Medium_Balloon = mean(Damping_Ratio_Medium_Balloon);

for n = 1:length(Volt_Long_Balloon)-1
    Period_Long = Time_Long_Balloon(Volt_Locations_Balloon_Long(n+1))-Time_Long_Balloon(Volt_Locations_Balloon_Long(n));
    Damping_Ratio_Long_Balloon(n) = 2*pi/Period_Long;
end
Damp_Mean_Long_Balloon = mean(Damping_Ratio_Long_Balloon);



Natural_Freq_Small_Valve = Damp_Mean_Small/sqrt(1-Ratio_Valve_Small^2);
Natural_Freq_Medium_Valve = Damp_Mean_Medium/sqrt(1-Ratio_Valve_Medium^2);
Natural_Freq_Long_Valve = Damp_Mean_Long/sqrt(1-Ratio_Valve_Long^2);

Natural_Freq_Small_Balloon = Damp_Mean_Small_Balloon/sqrt(1-Ratio_Balloon_Small^2);
Natural_Freq_Medium_Balloon = Damp_Mean_Medium_Balloon/sqrt(1-Ratio_Balloon_Medium^2);
Natural_Freq_Long_Balloon = Damp_Mean_Long_Balloon/sqrt(1-Ratio_Balloon_Long^2);


% Ratio_Valve_Small = mean(Damp_Small_Valve);
% Ratio_Valve_Medium = mean(Damp_Medium_Valve);
% Ratio_Valve_Long = mean(Damp_Long_Valve);
% Ratio_Balloon_Small = mean(Damp_Small_Balloon);
% Ratio_Balloon_Medium = mean(Damp_Medium_Balloon);
% Ratio_Balloon_Long = mean(Damp_Long_Balloon);

KA = 2.46;
e = exp((-Ratio_Valve_Small*Natural_Freq_Small_Valve.*Time_Small_Flood));
a = Ratio_Valve_Small/(sqrt(1-Ratio_Valve_Small^2));
b = Natural_Freq_Small_Valve.*Time_Small_Flood*sqrt(1-Ratio_Valve_Small^2);
Volt_Predict_Valve_Small = KA-KA.*e.*a.*sin(b)+cos(b);

figure(13)
plot(Time_Small_Flood,Volt_Predict_Valve_Small)


KA = 2.46;
e = exp((-Ratio_Valve_Medium*Natural_Freq_Medium_Valve.*Time_Medium_Flood));
a = Ratio_Valve_Medium/(sqrt(1-Ratio_Valve_Medium^2));
b = Natural_Freq_Medium_Valve.*Time_Medium_Flood*sqrt(1-Ratio_Valve_Medium^2);
Volt_Predict_Valve_Medium = KA-KA.*e.*a.*sin(b)+cos(b);

figure(14)
plot(Time_Medium_Flood,Volt_Predict_Valve_Medium)


KA = 2.46;
e = exp((-Ratio_Valve_Long*Natural_Freq_Long_Valve.*Time_Long_Flood));
a = Ratio_Valve_Long/(sqrt(1-Ratio_Valve_Long^2));
b = Natural_Freq_Long_Valve.*Time_Long_Flood*sqrt(1-Ratio_Valve_Long^2);
Volt_Predict_Valve_Long = KA-KA.*e.*a.*sin(b)+cos(b);

figure(15)
plot(Time_Long_Flood,Volt_Predict_Valve_Long)




KA = -.012;
e = exp((-Ratio_Balloon_Small*Natural_Freq_Small_Balloon.*Time_Small_Balloon));
a = Ratio_Balloon_Small/(sqrt(1-Ratio_Balloon_Small.^2));
b = Natural_Freq_Small_Balloon.*Time_Small_Balloon*sqrt(1-Ratio_Balloon_Small^2);
Volt_Predict_Balloon_Small = KA-KA.*e.*a.*sin(b)+cos(b);

figure(16)
plot(Time_Small_Balloon,Volt_Predict_Balloon_Small,Time_Small_Balloon,Voltage_Small_Balloon)


KA = -.01;
e = exp((-Ratio_Balloon_Medium*Natural_Freq_Medium_Balloon.*Time_Medium_Balloon));
a = Ratio_Balloon_Medium/(sqrt(1-Ratio_Balloon_Medium^2));
b = Natural_Freq_Medium_Balloon.*Time_Medium_Balloon*sqrt(1-Ratio_Balloon_Medium^2);
Volt_Predict_Balloon_Medium = KA-KA.*e.*a.*sin(b)+cos(b);

figure(17)
plot(Time_Medium_Balloon,Volt_Predict_Balloon_Medium)


KA = -.01;
e = exp((-Ratio_Balloon_Long*Natural_Freq_Long_Balloon.*Time_Long_Balloon));
a = Ratio_Balloon_Long/(sqrt(1-Ratio_Balloon_Long^2));
b = Natural_Freq_Long_Balloon.*Time_Long_Balloon*sqrt(1-Ratio_Balloon_Long^2);
Volt_Predict_Balloon_Long = KA-KA.*e.*a.*sin(b)+cos(b);

figure(18)
plot(Time_Long_Balloon,Volt_Predict_Balloon_Long)%,Time_Long_Balloon,Voltage_Long_Balloon)


