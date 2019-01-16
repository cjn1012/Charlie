clear all
close all

% This was written as a script for the ease of convience in grading. 
% The parameters below assume a tapered (small diameter at the bottom) cup with a hole in the center of the bottom
% with a circular cross section (circle hole)

g = 9.81; % m/s^2 , gravitational constant

% Parameters Needed assuming a cylindrical cup
WaterHeightInitial = 1;  % Meters
CupDiameterBottom = .5;  % Meters
CupDiameterTop = 1;      % Meters
HoleDiameter = .1;       % Meters

% Calculation to get the time
Area_Cup = pi * ((CupDiameterTop+CupDiameterBottom)/4)^2; % Area of cup with a taper accounting for converting to radii
Desired_Height_Half = WaterHeightInitial * 0.5;
Desired_Height_Quarter = WaterHeightInitial * 0.25;
Desired_Height_Eighth = WaterHeightInitial * 0.125;
Desired_Height_Sixteenth = WaterHeightInitial * 0.0625;
Area_Hole = pi * (HoleDiameter/2)^2;

% Equation derived from bernoulli and calculating time to reach a certain fraction of where it started
Time_to_Half = 2 * Area_Cup * (sqrt(WaterHeightInitial) - sqrt(Desired_Height_Half)) / (Area_Hole * sqrt(2*g));
Time_to_Quarter = 2 * Area_Cup * (sqrt(WaterHeightInitial) - sqrt(Desired_Height_Quarter)) / (Area_Hole * sqrt(2*g));
Time_to_Eighth = 2 * Area_Cup * (sqrt(WaterHeightInitial) - sqrt(Desired_Height_Eighth)) / (Area_Hole * sqrt(2*g));
Time_to_Sixteenth = 2 * Area_Cup * (sqrt(WaterHeightInitial) - sqrt(Desired_Height_Sixteenth)) / (Area_Hole * sqrt(2*g));

% Creating a table to display the data nicely 
Desired_Fraction      = [0.5;0.25;0.125;0.0625];
Initial_Water_Height  = [1;1;1;1];
Cup_Bottom_Diameter   = [0.5,0.5,0.5,0.5]';
Cup_Top_Diameter      = [1;1;1;1];
Hole_Diameter         = [0.1;0.1;0.1;0.1];
Time_to_Level         = [Time_to_Half,Time_to_Quarter,Time_to_Eighth,Time_to_Sixteenth]';
String  = 'The values below do make sense because the time increases as it needs to get more water out';
String2 = 'You also see the difference between the times to descrease as it needs to get less water out';
String3 = 'But it also is not directly linear, as it doesnt percicely get cut in half because the flow';
String4 = 'of water out is decrease as the height of the water decreases through time';


% Create a table, T, in Data_Test and then some sentences from above in 'thoughts'
T = table(Desired_Fraction,Initial_Water_Height,Cup_Bottom_Diameter,Cup_Top_Diameter,Hole_Diameter,Time_to_Level)
fid = fopen('Thoughts.txt','wt'); % Open text file for strings inserted
fprintf(fid,'%s\n%s\n%s\n%s', String, String2, String3,String4); % inserting string from above
fclose(fid);
writetable(T,'Data_Test.txt')



% Graphing the water level with time
X  = [0,0.25,0.45];
X2 = [0.55,0.75,1];
Y  = [1,0,0]
Y2 = [0,0,1]
hold on
plot(X,Y,'r')
plot(X2,Y2,'r')
text(0.35,1,'Full - t = 0','Color','b')
text(0.35,0.5,'Half - t = 7.44','Color','b')
text(0.35,0.25,'Quarter - t = 12.70','Color','b')
text(0.35,0.125,'Eighth - t = 16.41','Color','b')
text(0.35,0.0625,'Sixteenth - t = 19.05','Color','b')
xlim([-0.2,1.2])
ylim([-.2,1.2])
xlabel('Distance (m)')
ylabel('Distance (m)')











