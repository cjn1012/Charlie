% This was written as a script for the ease of convience in grading. 
% The parameters below assume a cylindrical cup with a hole in the center of the bottom
% with a circular cross section (circle hole)


g = 9.81; % m/s^2 , gravitational constant

% Parameters Needed assuming a cylindrical cup
WaterHeightInitial = 1;  % Meters
CupDiameter = 1;         % Meters
HoleDiameter = .1;       % Meters

% Calculation to get the time
Area_Cup = pi * (CupDiameter/2)^2; % Area of cup with a taper accounting for converting to radii
Desired_Height_Half = WaterHeightInitial * 0.5;
Area_Hole = pi * (HoleDiameter/2)^2;

% Equation derived from bernoulli and calculating time to reach a certain fraction of where it started
Time_to_Half = 2 * Area_Cup * (sqrt(WaterHeightInitial) - sqrt(Desired_Height_Half)) / (Area_Hole * sqrt(2*g));

% Creating a table to display the data nicely 
Desired_Fraction      = [0.5];
Initial_Water_Height  = [1];
Cup_Diameter          = [0.5];
Hole_Diameter         = [0.1];
Time_to_Level         = [Time_to_Half];

T = table(Desired_Fraction,Initial_Water_Height,Cup_Diameter,Hole_Diameter,Time_to_Level)

writetable(T,'Data_Test.txt')










