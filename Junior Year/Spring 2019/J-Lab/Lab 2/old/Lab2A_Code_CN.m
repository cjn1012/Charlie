clear all
close all

%% SECTION 1 - STATIC CALIBRATION



%%%%%%%%%%%   %%%%%%%%%%%  %%%%%%%%%%%  %%%%%%%%%%%%%  %%%%%%%%%%%%%  %%%%%%%%%%%  %        %       %%%%%%%%%%%% 
%             %            %                  %              %        %         %  % %      %            %                           
%             %            %                  %              %        %         %  %  %     %            %                                  
%             %            %                  %              %        %         %  %   %    %            %                          
%%%%%%%%%%%   %%%%%%%%%%%  %                  %              %        %         %  %    %   %            %                         
          %   %            %                  %              %        %         %  %     %  %            %                         
          %   %            %                  %              %        %         %  %      % %            %                 
          %   %            %                  %              %        %         %  %       %%            %                     
%%%%%%%%%%%   %%%%%%%%%%%  %%%%%%%%%%%        %        %%%%%%%%%%%%%  %%%%%%%%%%%  %        %       %%%%%%%%%%%%




%% Data Input - Hand Recorded

Zero_Temp     = [0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ]; % Degrees Celcius
Zero_Volt     = [.0347,.0383,.0362,.0314,.0366,.0372,.0337,.0318,.0339,.0355]; % Volt
Zero_Volt_Avg = mean(Zero_Volt); % Volt

Bath_Temp = [0            ,20   ,30   ,40   ,50   ,60   ,70   ,80   ,90   ,100  ]; % Degrees Celcius
Ther_Resi = [27200        ,12190,8160 ,5650 ,3940 ,2890 ,2160 ,1720 ,1340 ,1020 ]; % Ohms
Amp_Volt  = [Zero_Volt_Avg,.2317,.3285,.4321,.5387,.6416,.7441,.8171,.9484,1.020]; % Volt


% Done


%% Part 1 - Done within Word

% Done

%% Part 3 - Calibration Linear Curve - Using Resistance to Calculate Temperature of the *Thermistor*


Beta = 3343; % (K)   - Calculated on Document
R_0  = 9736; % (Ohm) - Calculated on Document
T_0  = 298.15 ; % (K)   - Given as Room Temperature Water

Temp_Thermistor_Actual = zeros(10,1); % Creating array of zeros

for i = 1:length(Bath_Temp) % Calculating the Temperature of the Thermistor from the Resistance Data
    Temp_Thermistor_Actual(i) =  ( ( (1/Beta) * log((Ther_Resi(i))/R_0) + (1/T_0) )^-1 ) -273.15; 
end

Fit_Voltage = polyfit(Temp_Thermistor_Actual,Amp_Volt',1); % Calculates the slope and y-intercept to calculate fitted data of Voltage from actual temperatures
Fit_1_Voltage = Fit_Voltage(1).*Temp_Thermistor_Actual + Fit_Voltage(2); % Calculates the 10 voltages that use the paramaters calculated above

% Plotting of the raw Bath temperature and voltage measurement with the best fit line
figure(1)
plot(Temp_Thermistor_Actual,Fit_1_Voltage,'Color','b')
hold on
plot(Bath_Temp,Amp_Volt,'o','Color','r')
xlabel('Bath Temperature (°C)')
ylabel('Voltage (V)')
xlim([-5,105])
ylim([0,1.1])
text(50,0.4,'Voltage(Volts) = 0.01*T(°C) + 0.035')
title('Calibration Curve for the Thermocouple')
legend('Best Fit', 'Data Points', 'Location', 'Northwest')


% Done

%% Part 4 - Converting Thermocouple Voltage Measurements to Temperature. Compare with Thermistor Measurement

% Calculating the Thermocouple Temperatures in the Water baths from the calibration curve found in Part 3
N = length(Amp_Volt); 
Temp_Thermocouple_Actual = zeros(N,1);
for i = 1:N
    Temp_Thermocouple_Actual(i) = (Amp_Volt(i)-Fit_Voltage(2))/Fit_Voltage(1); % Calculates Temperature measurements from the voltages recorded
end

% Calculates the parameters for a fit line comparing the Thermocouple and Thermistor Temperatures
Best_Fit_2 = polyfit(Temp_Thermocouple_Actual,Temp_Thermistor_Actual,1); 

% Calculating the Values of the fit line comparing the Thermocouple to the Thermistor
Fit_2 = Best_Fit_2(1).*Temp_Thermistor_Actual + Best_Fit_2(2); %2


nu = length(Temp_Thermistor_Actual)-(2); % Calculating nu for statistic calculation
t_nup =tinv(.95, nu); % Calculating the t value for measuring confidence for fit data and measured
Temp_Zero_Mean = mean(Temp_Thermistor_Actual); % Calculating the mean of Thermistor
sumbot = sum((Temp_Thermistor_Actual-Temp_Zero_Mean).^2); % Denominator of the confidence calculation
syx = (sum((Temp_Thermocouple_Actual-Fit_2).^2)/nu).^.5; % Numerator for conidence calculation

% Calculating the confidence interval for the fit
confit = t_nup*syx*(1/length(Temp_Thermistor_Actual)+(Temp_Thermistor_Actual-Temp_Zero_Mean).^2/ sumbot).^.5; 
% Calculating the confidence interval for the measured
conmeas = t_nup*syx*(1+1/length(Temp_Thermistor_Actual)+(Temp_Thermistor_Actual-Temp_Zero_Mean).^2/sumbot).^.5;

% Plotting
figure(2)
plot(Temp_Thermistor_Actual,Temp_Thermocouple_Actual,'o')
hold on
plot(Temp_Thermistor_Actual,Fit_2)
plot(Temp_Thermistor_Actual, Fit_2+confit, Temp_Thermistor_Actual, Fit_2-confit, Temp_Thermistor_Actual, Fit_2+conmeas, Temp_Thermistor_Actual, Fit_2-conmeas)
xlabel('Thermistor Temperature Reading (^oC)') 
ylabel('Thermocouple Temperature Reading (^oC)') 
text(60,30,['t_9_5_%_,_8=' num2str(t_nup)])
title('Temperature Readings with Confidence Intervals of Fit and Measured') 
legend('Data Points', 'Best Fit', 'Upper Confidence Bound of Fit', 'Lower Confidence Bound of Fit', 'Upper Confidence Bound of Measurement', 'Lower Confidence Bound of Measurement',  'Location', 'northwest') 
xlim([-10,105])
ylim([-5,105])
grid minor


% Done


%% Part 2 - Repetitive Ice Water Measurements - Used for Word Document Submission

N = length(Zero_Temp);
Zero_Temp_Actual = zeros(N,1);
for i = 1:N
    Zero_Temp_Actual(i) = (Zero_Volt(i)-Fit_Voltage(2))/Fit_Voltage(1); % Calculates Temperature measurements from the voltages recorded
end

% Mean Calculation of Temperatures
Zero_Temp_Mean = mean(Zero_Temp_Actual);

for i = 1:length(Zero_Temp)
    Std_Dev = sqrt(((Zero_Temp_Actual(i)-Zero_Temp_Mean)^2)/(N-1)); 
end

% From Tables
t_vP = 2.262; % 95% Confidence, 10 Data Points

% Calculations
Std_Dev_Mean = Std_Dev / N;
Population_Mean_Upper = Zero_Temp_Mean + t_vP*Std_Dev_Mean;
Population_Mean_Lower = Zero_Temp_Mean - t_vP*Std_Dev_Mean;


% Done


%% Part 5 - Ice Bath Measurement Comparison to the fit and confidence interval of fit and measurements

% Plotting
figure(3)
plot(0,Zero_Temp_Actual,'o') % Data point of Calculated temperature and 0 degrees C
hold on
% Plotting confidence intervals from Part 4 on graph to compare
plot(Temp_Thermistor_Actual, Fit_2+confit, Temp_Thermistor_Actual, Fit_2-confit, Temp_Thermistor_Actual, Fit_2+conmeas, Temp_Thermistor_Actual, Fit_2-conmeas)
ylim([-3,3])
xlim([-.1,.5])
xlabel('Temperature (^oC)')
ylabel('Temperature (^oC)')
legend('Ice Bath Data Points', 'Upper Confidence Bound of Measurement', 'Lower Confidence Bound of Measurement', 'Upper Confidence Bound of Fit', 'Lower Confidence Bound of Fit',  'Location', 'northeast')
title('Ice Bath Temperature Calculation with Confidence Intervals')

% Done

%% SECTION 2 - DYNAMIC CALIBRATION


%%%%%%%%%%%   %%%%%%%%%%%  %%%%%%%%%%%  %%%%%%%%%%%%%  %%%%%%%%%%%%%  %%%%%%%%%%%  %        %       %%%%%%%%%%%%%%%%%%%%%  
%             %            %                  %              %        %         %  % %      %            %        %                    
%             %            %                  %              %        %         %  %  %     %            %        %                          
%             %            %                  %              %        %         %  %   %    %            %        %                  
%%%%%%%%%%%   %%%%%%%%%%%  %                  %              %        %         %  %    %   %            %        %                 
          %   %            %                  %              %        %         %  %     %  %            %        %                 
          %   %            %                  %              %        %         %  %      % %            %        %         
          %   %            %                  %              %        %         %  %       %%            %        %             
%%%%%%%%%%%   %%%%%%%%%%%  %%%%%%%%%%%        %        %%%%%%%%%%%%%  %%%%%%%%%%%  %        %       %%%%%%%%%%%%%%%%%%%%%

%% Part 1 

% BIB = Bare Ice Boil
% BBI = Bare Boil Ice
% BIA = Bare Ice Air
% BIW = Bare Ice Water
% AIB = Aluminum Ice Boil
% ABI = Aluminum Boil Ice
% SIB = Steel Ice Boil
% SBI = Steel Boil Ice

% Bare Thermocoupler, N = 50,000
BIB_Time = xlsread('Part2Data.xlsx','bareiceboil' ,'A9:A50008'); % Second
BIB_Volt = xlsread('Part2Data.xlsx','bareiceboil' ,'B9:B50008'); % Volt
BBI_Time = xlsread('Part2Data.xlsx','bareboilice' ,'A9:A50008');
BBI_Volt = xlsread('Part2Data.xlsx','bareboilice' ,'B9:B50008');

% Bare Thermocoupler, N = 12,000
BIA_Time = xlsread('Part2Data.xlsx','bareiceair'  ,'A9:A50008'); % Second
BIA_Volt = xlsread('Part2Data.xlsx','bareiceair'  ,'B9:B50008'); % Volt
BIW_Time = xlsread('Part2Data.xlsx','bareicewater','A9:A50008');
BIW_Volt = xlsread('Part2Data.xlsx','bareicewater','B9:B50008');


% Aluminum and Steel Insulated Thermocoupler, N = 5,000
AIB_Time = xlsread('Part2Data.xlsx','aliceboil'   ,'A9:A5008'); % Second
AIB_Volt = xlsread('Part2Data.xlsx','aliceboil'   ,'B9:B5008'); % Volt
ABI_Time = xlsread('Part2Data.xlsx','alboilice'   ,'A9:A5008');
ABI_Volt = xlsread('Part2Data.xlsx','alboilice'   ,'B9:B5008');
SIB_Time = xlsread('Part2Data.xlsx','steeliceboil','A9:A5008');
SIB_Volt = xlsread('Part2Data.xlsx','steeliceboil','B9:B5008');
SBI_Time = xlsread('Part2Data.xlsx','steelboilice','A9:A5008');
SBI_Volt = xlsread('Part2Data.xlsx','steelboilice','B9:B5008');


% Done



% Using the previous calibration curve to convert all voltage measurements to temperature

BIB_Temp = (BIB_Volt - Fit_Voltage(2))./Fit_Voltage(1);
BBI_Temp = (BBI_Volt - Fit_Voltage(2))./Fit_Voltage(1); 
BIA_Temp = (BIA_Volt - Fit_Voltage(2))./Fit_Voltage(1);
BIW_Temp = (BIW_Volt - Fit_Voltage(2))./Fit_Voltage(1);
AIB_Temp = (AIB_Volt - Fit_Voltage(2))./Fit_Voltage(1);
ABI_Temp = (ABI_Volt - Fit_Voltage(2))./Fit_Voltage(1);
SIB_Temp = (SIB_Volt - Fit_Voltage(2))./Fit_Voltage(1);
SBI_Temp = (SBI_Volt - Fit_Voltage(2))./Fit_Voltage(1);


% Smoothing data

Span   = 200;
Window = ones(Span,1)/Span;
BIB_Temp_Smooth = conv(BIB_Temp,Window,'same');
BIB_Temp_Smooth = BIB_Temp_Smooth(1000:end-1000);
BIB_Time = BIB_Time(1000:end-1000);

Span   = 200;
Window = ones(Span,1)/Span;
BBI_Temp_Smooth = conv(BBI_Temp,Window,'same');
BBI_Temp_Smooth = BBI_Temp_Smooth(1000:end-1000);
BBI_Time = BBI_Time(1000:end-1000);

Span   = 100;
Window = ones(Span,1)/Span;
BIA_Temp_Smooth = conv(BIA_Temp,Window,'same');
BIA_Temp_Smooth = BIA_Temp_Smooth(250:end-250);
BIA_Time = BIA_Time(250:end-250);

Span   = 100;
Window = ones(Span,1)/Span;
BIW_Temp_Smooth = conv(BIW_Temp,Window,'same');
BIW_Temp_Smooth = BIW_Temp_Smooth(250:end-250);
BIW_Time = BIW_Time(250:end-250);

Span   = 50;
Window = ones(Span,1)/Span;
AIB_Temp_Smooth = conv(AIB_Temp,Window,'same');
AIB_Temp_Smooth = AIB_Temp_Smooth(100:end-100);
AIB_Time = AIB_Time(100:end-100);

Span   = 50;
Window = ones(Span,1)/Span;
ABI_Temp_Smooth = conv(ABI_Temp,Window,'same');
ABI_Temp_Smooth = ABI_Temp_Smooth(100:end-100);
ABI_Time = ABI_Time(100:end-100);

Span   = 50;
Window = ones(Span,1)/Span;
SIB_Temp_Smooth = conv(SIB_Temp,Window,'same');
SIB_Temp_Smooth = SIB_Temp_Smooth(100:end-100);
SIB_Time = SIB_Time(100:end-100);

Span   = 50;
Window = ones(Span,1)/Span;
SBI_Temp_Smooth = conv(SBI_Temp,Window,'same');
SBI_Temp_Smooth = SBI_Temp_Smooth(100:end-100);
SBI_Time = SBI_Time(100:end-100);


% Using a sliding, first order polynomial to determine time of maximum slope

%BIB
Window = 15;
MaxSlope = 0;
for i = Window+1:length(BIB_Time)-Window-1
    p = polyfit(BIB_Time(i-Window:i+Window),BIB_Temp_Smooth(i-Window:i+Window),1);
    if (abs(p(1))>MaxSlope)
        MaxSlope = abs(p(1));
        startfitBIB = i;
    end
end
BIB_Time = BIB_Time - BIB_Time(startfitBIB);


%BBI
Window = 15;
MaxSlope = 0;
for i = Window+1:length(BBI_Time)-Window-1
    p = polyfit(BBI_Time(i-Window:i+Window),BBI_Temp_Smooth(i-Window:i+Window),1);
    if (abs(p(1))>MaxSlope)
        MaxSlope = abs(p(1));
        startfitBBI = i;
    end
end
BBI_Time = BBI_Time - BBI_Time(startfitBBI);


%BIA
Window = 10;
MaxSlope = 0;
for i = Window+1:length(BIA_Time)-Window-1
    p = polyfit(BIA_Time(i-Window:i+Window),BIA_Temp_Smooth(i-Window:i+Window),1);
    if (abs(p(1))>MaxSlope)
        MaxSlope = abs(p(1));
        startfitBIA = i;
    end
end
BIA_Time = BIA_Time - BIA_Time(startfitBIA);

%BIW
Window = 10;
MaxSlope = 0;
for i = Window+1:length(BIW_Time)-Window-1
    p = polyfit(BIW_Time(i-Window:i+Window),BIW_Temp_Smooth(i-Window:i+Window),1);
    if (abs(p(1))>MaxSlope)
        MaxSlope = abs(p(1));
        startfitBIW = i;
    end
end
BIW_Time = BIW_Time - BIW_Time(startfitBIW);

%AIB
Window = 5;
MaxSlope = 0;
for i = Window+1:length(AIB_Time)-Window-1
    p = polyfit(AIB_Time(i-Window:i+Window),AIB_Temp_Smooth(i-Window:i+Window),1);
    if (abs(p(1))>MaxSlope)
        MaxSlope = abs(p(1));
        startfitAIB = i;
    end
end
AIB_Time = AIB_Time - AIB_Time(startfitAIB);


%ABI
Window = 5;
MaxSlope = 0;
for i = Window+1:length(ABI_Time)-Window-1
    p = polyfit(ABI_Time(i-Window:i+Window),ABI_Temp_Smooth(i-Window:i+Window),1);
    if (abs(p(1))>MaxSlope)
        MaxSlope = abs(p(1));
        startfitABI = i;
    end
end
ABI_Time = ABI_Time - ABI_Time(startfitABI);


%SIB
Window = 5;
MaxSlope = 0;
for i = Window+1:length(SIB_Time)-Window-1
    p = polyfit(SIB_Time(i-Window:i+Window),SIB_Temp_Smooth(i-Window:i+Window),1);
    if (abs(p(1))>MaxSlope)
        MaxSlope = abs(p(1));
        startfitSIB = i;
    end
end
SIB_Time = SIB_Time - SIB_Time(startfitSIB);

%SBI
Window = 5;
MaxSlope = 0;
for i = Window+1:length(SBI_Time)-Window-1
    p = polyfit(SBI_Time(i-Window:i+Window),SBI_Temp_Smooth(i-Window:i+Window),1);
    if (abs(p(1))>MaxSlope)
        MaxSlope = abs(p(1));
        startfitSBI = i;
    end
end
SBI_Time = SBI_Time - SBI_Time(startfitSBI);



% Plotting

figure(4)
plot(BIB_Time,BIB_Temp_Smooth)
hold on
plot(BIB_Time(startfitBIB),BIB_Temp_Smooth(startfitBIB),'o')
title('Bare Thermocouple Ice Bath to Boiling')
xlabel('Time (s)')
ylabel('Temperature (°C)')
ylim([-10,110])
xlim([-2,2])
text(BIB_Time(startfitBIB)+.1,BIB_Temp_Smooth(startfitBIB),strcat('T_{initial} = ' ,num2str(BIB_Temp_Smooth(startfitBIB)),'°C'))


figure(5)
plot(BBI_Time,BBI_Temp_Smooth)
hold on
plot(BBI_Time(startfitBBI),BBI_Temp_Smooth(startfitBBI),'o')
title('Bare Thermocouple Boiling to Ice Bath')
xlabel('Time (s)')
ylabel('Temperature (°C)')
ylim([-10,110])
xlim([-2,2])
text(BBI_Time(startfitBBI)+.1,BBI_Temp_Smooth(startfitBBI),strcat('T_{initial} = ' ,num2str(BBI_Temp_Smooth(startfitBBI)),'°C'))


figure(6)
plot(AIB_Time,AIB_Temp_Smooth)
hold on
plot(AIB_Time(startfitAIB),AIB_Temp_Smooth(startfitAIB),'o')
title('Aluminum-Insulated Thermocouple Ice Bath to Boiling')
xlabel('Time (s)')
ylabel('Temperature (°C)')
ylim([-10,110])
xlim([-10,30])
text(AIB_Time(startfitAIB)+1,AIB_Temp_Smooth(startfitAIB),strcat('T_{initial} = ' ,num2str(AIB_Temp_Smooth(startfitAIB)),'°C'))


figure(7)
plot(ABI_Time,ABI_Temp_Smooth)
hold on
plot(ABI_Time(startfitABI),ABI_Temp_Smooth(startfitABI),'o')
title('Aluminum-Insulated Thermocouple Boiling to Ice Bath')
xlabel('Time (s)')
ylabel('Temperature (°C)')
ylim([-10,110])
xlim([-10,40])
text(ABI_Time(startfitABI)+1,ABI_Temp_Smooth(startfitABI),strcat('T_{initial} = ' ,num2str(ABI_Temp_Smooth(startfitABI)),'°C'))


figure(8)
plot(SIB_Time,SIB_Temp_Smooth)
hold on
plot(SIB_Time(startfitSIB),SIB_Temp_Smooth(startfitSIB),'o')
title('Steel-Insulated Thermocouple Ice Bath to Boiling')
xlabel('Time (s)')
ylabel('Temperature (°C)')
ylim([-10,110])
xlim([-5,20])
text(SIB_Time(startfitSIB)+1,SIB_Temp_Smooth(startfitSIB),strcat('T_{initial} = ' ,num2str(SIB_Temp_Smooth(startfitSIB)),'°C'))


figure(9)
plot(SBI_Time,SBI_Temp_Smooth)
hold on
plot(SBI_Time(startfitSBI),SBI_Temp_Smooth(startfitSBI),'o')
title('Steel-Insulated Thermocouple Boiling to Ice Bath')
xlabel('Time (s)')
ylabel('Temperature (°C)')
ylim([-10,110])
xlim([-5,40])
text(SBI_Time(startfitSBI)+1,SBI_Temp_Smooth(startfitSBI),strcat('T_{initial} = ' ,num2str(SBI_Temp_Smooth(startfitSBI)),'°C'))


%% Part 2

% BIB = Bare Ice Boil
% BBI = Bare Boil Ice
% BIA = Bare Ice Air
% BIW = Bare Ice Water
% AIB = Aluminum Ice Boil
% ABI = Aluminum Boil Ice
% SIB = Steel Ice Boil
% SBI = Steel Boil Ice

% % Calculate Final Temperature Average for each Curve
BIB_Temp_Final = mean(BIB_Temp_Smooth(end-1000:end));
BBI_Temp_Final = mean(BBI_Temp_Smooth(end-1000:end));
AIB_Temp_Final = mean(AIB_Temp_Smooth(end-100:end));
ABI_Temp_Final = mean(ABI_Temp_Smooth(end-100:end));
SIB_Temp_Final = mean(SIB_Temp_Smooth(end-100:end));
SBI_Temp_Final = mean(SBI_Temp_Smooth(end-100:end));







