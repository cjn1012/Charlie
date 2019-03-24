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

a = 343;
l = .17;
V = l*pi*.0042^2;
Vt = 6.5*10^-8;
Damp_Natural_Freq_Small_Flood = a/(l*sqrt(0.5+(V/Vt)));

figure(1)

subplot(2,1,1)
plot (Time_Small_Flood,Voltage_Small_Flood);
title('Sine vs. Triangle Wave')
xlabel ('Time (sec)')
ylabel ('Volts (V)')

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

a = 343;
l = .795;
V = l*pi*.0042^2;
Vt = 6.5*10^-8;
Damp_Natural_Freq_Medium_Flood = a/(l*sqrt(0.5+(V/Vt)));

figure(2)

subplot(2,1,1)
plot (Time_Medium_Flood,Voltage_Medium_Flood);
title('Sine vs. Triangle Wave')
xlabel ('Time (sec)')
ylabel ('Volts (V)')

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

a = 343;
l = 1.05;
V = l*pi*.0042^2;
Vt = 6.5*10^-8;
Damp_Natural_Freq_Long_Flood = a/(l*sqrt(0.5+(V/Vt)));

figure(3)

subplot(2,1,1)
plot (Time_Long_Flood,Voltage_Long_Flood);
title('Sine vs. Triangle Wave')
xlabel ('Time (sec)')
ylabel ('Volts (V)')

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

a = 343;
l = .17;
V = l*pi*.0042^2;
Vt = 6.5*10^-8;
Damp_Natural_Freq_Small_Balloon = a/(l*sqrt(0.5+(V/Vt)));

figure(4)

subplot(2,1,1)
plot (Time_Small_Balloon,Voltage_Small_Balloon);
title('Sine vs. Triangle Wave')
xlabel ('Time (sec)')
ylabel ('Volts (V)')

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

a = 343;
l = .795;
V = l*pi*.0042^2;
Vt = 6.5*10^-8;
Damp_Natural_Freq_Medium_Balloon = a/(l*sqrt(0.5+(V/Vt)));

figure(5)

subplot(2,1,1)
plot (Time_Medium_Balloon,Voltage_Medium_Balloon);
title('Sine vs. Triangle Wave')
xlabel ('Time (sec)')
ylabel ('Volts (V)')

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

a = 343;
l = 1.05;
V = l*pi*.0042^2;
Vt = 6.5*10^-8;
Damp_Natural_Freq_Long_Balloon = a/(l*sqrt(0.5+(V/Vt)));

figure(6)

subplot(2,1,1)
plot (Time_Long_Balloon,Voltage_Small_Balloon);
title('Sine vs. Triangle Wave')
xlabel ('Time (sec)')
ylabel ('Volts (V)')

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







