% Get a randomized initial conditions for the trajectoryNumber defined
initialConditions = GetAllInitialConditions(fromInitialCondition, toInitialCondition, trajectoryNumber);

allTrajectories = GetAllTrajectories(initialConditions, t);

%%

% Plot 3D view of the potential landscape
figure;

% Get the probability distributions
positionProbabilities = GetPositionProbabilities(fromInitialCondition, toInitialCondition, index, splitNumber, allTrajectories);

% Normalized Probability
normalized_Prob=positionProbabilities./(sum(sum(positionProbabilities)));

% U = - ln P
Neg_lnP=-log(normalized_Prob); 

[x y] = meshgrid(linspace(fromInitialCondition(index(1)), toInitialCondition(index(1)), splitNumber), ...
    linspace(fromInitialCondition(index(2)), toInitialCondition(index(2)), splitNumber));

% plot the 3D view of the potential landscape
surf(x,y,Neg_lnP);

% label the axes
xlabel(variableNames(index(2)));
ylabel(variableNames(index(1)));
zlabel('U')

% set the x and y axes
axis([range_min range_max range_min range_max]);

% 3D view
% for Figure 6B view(-14,61);
view(az_angle,el_angle); 
box on;

% set the color bar 
colorbar;
colormap(jet);

%%

% Plot Top view of the potential landscape
figure;

surf(x,y,Neg_lnP);
xlabel(variableNames(index(2)));
ylabel(variableNames(index(1)));
zlabel('U')

axis([range_min range_max range_min range_max]);

% Top view
view(0,90);
box on;
colorbar;
% for MATLAB R2017b you need to add below colormap with jet for color from
% red (large value) to blue (small value)
colormap(jet);

% show how long it takes to calculate and plot the landscape
disp('Calculation and plotting finished.')
duration=toc/(60*60);
disp(['It took ' num2str(duration) ' hours'])

