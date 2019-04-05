textFileName = 'Data.txt';
% This script selects a random number from a text file with a list of numbers in separated in columns

% This script requires an input that is the file name of a text file in the working directory.
% This filename should be inserted as a string without the '.txt'.

numbers = csvread(textFileName);    % Pulls in desired text file as numbers array
randomPos = randi(length(numbers)); % Creates a random index between one and the length of numbers
randomNumber = numbers(randomPos)   % Assigns the answer, random number, to the numbers array indexed at randomPos

