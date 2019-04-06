% Engine Lab
clear all
close all

%% Reading in Data
% Reading in entire data set
Header=22;
RPM600  = importdata('600RPM.lvm' ,'\t',Header);
RPM700  = importdata('700RPM.lvm' ,'\t',Header);
RPM800  = importdata('800RPM.lvm' ,'\t',Header);
RPM900  = importdata('900RPM.lvm' ,'\t',Header);
RPM1000 = importdata('1000RPM.lvm','\t',Header);
RPM1100 = importdata('1100RPM.lvm','\t',Header);

% Separating entire data sets into vectors of Time and Voltage
RPM600_Volt  = RPM600.data(:,2); % Optical Sensor (V)
RPM700_Volt  = RPM700.data(:,2);
RPM800_Volt  = RPM800.data(:,2);
RPM900_Volt  = RPM900.data(:,2);
RPM1000_Volt = RPM1000.data(:,2);
RPM1100_Volt = RPM1100.data(:,2); 

RPM600_Pres  = RPM600.data(:,3); % p(V)
RPM700_Pres  = RPM700.data(:,3);
RPM800_Pres  = RPM800.data(:,3);
RPM900_Pres  = RPM900.data(:,3);
RPM1000_Pres = RPM1000.data(:,3);
RPM1100_Pres = RPM1100.data(:,3);


%% Contants for the Engine

vp_conv  = 0.0104; % V/psi
CR       = 8.5; % Compression Ratio
Rod_L    = .116; % Meters
Stroke   = .067; % Meters
Vol_Disp = 624/(1e6); % Cubic Meters
Vol_Clea = Vol_Disp/(CR-1); % Cubic Meters
R        = 2*Rod_L/Stroke; 

%% 600 RPM

RPM600P = (RPM600_Pres./vp_conv).*6894; % Pa
RPM600t = zeros(1,length(RPM600P));
















