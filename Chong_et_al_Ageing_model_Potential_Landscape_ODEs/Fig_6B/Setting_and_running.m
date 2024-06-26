%% Setting_and_running.m implements the algorithm [1] to calculate the potential landscape, plot the figures in 3D view and Top view, 
% and finally record the time it takes.
% 1. Zhang, X., Chong, K.H.  & Zheng, J. (2018). A Monte Carlo method for in silico modeling and visualization of
%?    Waddington's epigenetic landscape with intermediate details.


% Author: Ket Hing Chong
% School of Computer Science and Engineering
% Nanyang Technological University
% Singapore
% Email: kething@yahoo.com.sg
% Last revision: 7 May 2018

clear all;
close all;

addpath('../common_code'); % to access the common MATLAB files
% time in clock to record the start time
tic;

% define the variables name using the order as in ODEs equations
% For the example 1, there are 2 variables
variableNames = {'P53', 'Mdm2', 'Wip1', 'ATM', 'P21', 'PTEN', 'AKT', 'Myc', 'E2F', 'RB', 'CycE', 'CycD', 'ARF'};

range_max=5;
range_min=0;
% set the range of the space of interest (you can get maximum range from time course simulations)
initialRange = [zeros(1,13)+range_min;    % set the minimum range
               zeros(1,13)+range_max]; % set the maximum range

index = [1 4]; % the index of the two variables of interest to view the phase plane (1 for P53 and 2 for ATM)
           
% t define the numerical integration interval
t = 0:0.1:15; % set the max end time           
           
fromInitialCondition = min(initialRange);
toInitialCondition = max(initialRange);

% set the number of trajectories for time course simulations 
trajectoryNumber = 100000; 

% set the number of grid box size (100 x 100)
splitNumber = 100;

% 3D view setting
az_angle = -14;
el_angle = 61;
% Call the DrawLandscape function to draw landscape:
DrawLandscape;