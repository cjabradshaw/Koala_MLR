program PopulationSize_analyses

% the script evaluates koala population trends and the financial implications of different sterilisation strategies, 
% helping to identify the most cost-effective approach for population control.

% This script performs a multi-step analysis of koala population size and the costs of different sterilisation strategies over time. 
% Here's a brief summary of the key steps:
%%   Population Size Analysis (STEP 0): 
%       The script visualizes the female koala population over time, comparing the baseline scenario with Scenario 0. It uses shaded areas to represent the uncertainty between low and high confidence intervals.
%%   Baseline Scenario without Sterilisation (STEP 1):
%       The script simulates and plots the koala population's demographic evolution without any sterilisation (unmanaged scenario). It assesses how koala density compares to a threshold of 70 koalas per hectare and redistributes populations based on habitat suitability. The future population is evaluated, and areas where density exceeds the threshold are identified and saved for future analysis.
%%   Impact of Sterilisation on Population (STEP 2):
%       This step models the impact of sterilisation strategies (Scenario 1: female-only sterilisation and Scenario 2: female + daughter sterilisation) on the total koala population size. The script uses a range of sterilisation proportions (from 0 to 0.9) to predict future population changes under these scenarios, with reference lines showing population thresholds. It plots the outcomes of each scenario and highlights the sterilisation threshold where acceptable population control is achieved.
%%   Cost Calculation Based on Sterilisation Proportion (STEP 3):
%       The script calculates the cost of sterilisation for both scenarios (female-only and female + daughter) over time. It interpolates the cost data and visualizes it using color maps, where each map shows the cost as a function of sterilisation proportion and time. The maximum cost across both scenarios is used to define the color scale.
%%   Total Cost Over Time (STEP 4):
%       Finally, the script calculates the total cost of sterilisation over 25 years for each scenario by finding the sterilisation proportion closest to the desired value and summing the costs. The results provide a comparison between the total costs of the two strategies.
% author Frederik Saltre 13/09/2024

%%=========================================================================
%% STEP 0: DOUBLE CHECKING THAT DIFFERENT MODEL DO THE RIGHT THING
%%=========================================================================
% This setup provides a visualization of the total female population size over time, 
% showing both the baseline scenario and Scenario 0, with shading representing confidence intervals for uncertainty.

clear all, close all, clc,% Clear workspace, close all figures, and clear the command window

% Read baseline median values from CSV file
mat0 = readmatrix('Unmanaged_Scenario(Median).csv'); % Baseline scenario (median) values

% Read median, high CI, and low CI values for total female population over time for Scenario 0
fmed0 = readmatrix('Total_popFemale(Median)_OverTime_scenario#0.csv'); % Scenario 0 median values
fup0 = readmatrix('Total_popFemale(HighCI)_OverTime_Scenario#0.csv'); % Scenario 0 high confidence interval values
flo0 = readmatrix('Total_popFemale(LowCI)_OverTime_scenario#0.csv'); % Scenario 0 low confidence interval values

[r,c]=size(fmed0);% Get the size (number of rows and columns) of the median value data

% Create a new figure (figure 1) to plot the data
% Shade the area between two columns of 'mat0' representing uncertainty for the baseline scenario
% The columns represent the years and the range between the third and fourth column is filled with blue color
% Shade the area between the low and high CI values for Scenario 0 (green color)
figure(1),shade(mat0(:,1), mat0(:,3), mat0(:,1), mat0(:,4), 'FillType', [1, 2; 2, 1], 'FillColor', 'blue', 'FillAlpha', 0.3);
          hold on, shade(mat0(:,1), flo0(2,2:c), mat0(:,1), fup0(2,2:c), 'FillType', [1, 2; 2, 1], 'FillColor', 'green', 'FillAlpha', 0.3);
          % Plot the baseline median values (blue line) using the second column of 'mat0'
          % Plot the median values for the total female population in Scenario 0 (green line)
          hold on, plot(mat0(:,1), mat0(:,2),'-b','LineWidth',1.5),hold on, plot(mat0(:,1), fmed0(2,2:c),'-g','LineWidth',1.5),
          axis square; xlab('years', 13);ylab('total population (female) size (N)',13);% Set axis to square shape and label the x-axis and y-axis



%%=========================================================================
%% STEP 1: THE BASELINE = NO STERILISATION
%%=========================================================================
% This section plots the demographic evolution of the koala population if there is no fertility control
% Explanation of Key Sections:
%   Scenario Data Loading: mat0, mat1, and mat2 contain the baseline, low CI, and high CI data, respectively. The data is scaled by multiplying by 2.
%   Plotting Baseline Scenario: The shade function is used to fill the area between the low and high CI values for the baseline scenario in blue.
%   Koala Density Evaluation: The script evaluates koala densities based on habitat suitability and checks if they exceed a threshold of 70 koalas per km².
%                             It calculates the total number of koalas above the threshold in the present and the future.
%   Redistribution Based on Habitat Suitability: Population is redistributed based on the habitat suitability percentages. This is done both for present and future populations.
%   Threshold Check and Save Data: The script checks how many locations have koala densities exceeding 70 per hectare in both present and future scenarios. It saves this data for further analysis or plotting.
%   Percentage of Grid Cells Above Threshold: The script calculates the percentage of grid cells with koala populations above the threshold for both present and future scenarios.

clear all, close all, clc,% Clear workspace, close all figures, and clear the command window
cd '~/....'% Change directory to the folder where unmanaged scenario data is stored. To be modified as needed

% Load data for the unmanaged scenario (no intervention)
mat0 = readmatrix('Unmanaged_Scenario(Median).csv');mat0(:,2:4)=mat0(:,2:4).*2;% Load median values + Scale median values by 2
mat1 = readmatrix('Unmanaged_Scenario(LowCI).csv');mat1(:,2:4)=mat1(:,2:4).*2;% Load low confidence interval values + Scale low CI values by 2
mat2 = readmatrix('Unmanaged_Scenario(HighCI).csv');mat2(:,2:4)=mat2(:,2:4).*2;% Load high confidence interval values + Scale high CI values by 2

stats=[median(mat0(:,2:4),1);median(mat1(:,2:4),1);median(mat2(:,2:4),1)];% Calculate the median of the scaled values for the demographic statistics

% Plot the demographic evolution (shading between the low and high CI values)
figure(1),shade(mat0(:,1), mat1(:,3), mat0(:,1), mat2(:,4), 'FillType', [1, 2; 2, 1], 'FillColor', 'blue', 'FillAlpha', 0.3);


%% Evaluating how many koalas are and will be over the DEW conservation plan (0.7 koalas/ha)
% We work with the median values here. The density threshold is set to 70 koalas per km².
% The goal is to determine areas where koala density exceeds this threshold.
% M =>  [,1] = lon ; [,2]= lat ; [,3] = habitat_suitability ;	[,4] = mean_density ; [,5] = min_density ; [,6] = max_density
% we check for each column [4:6] = how much the densitu is over the treshold of 70 koala per km2

cd '~/....'% Change directory to the main folder for the koala model runs. To be modified as needed

M = readmatrix('Koala_density-SDM_outputs.csv');tr=70;% Load habitat suitability and koala density data from a CSV file + Density threshold (70 koalas per km²)

out=M(:,4:6)-tr;[r,c]=find(out<0);out2=out;% Subtract the threshold value from each column [4:6] (density values)
% Find locations where density is below the threshold and replace those values with NaN
for i=1:length(r);out2(r(i),c(i))=NaN;end;

% Calculate the total number of koalas above the threshold (omit NaN values)
stat70=sum(out2,1,"omitnan"); % This gives the total koala population above the threshold
% Example result: stat70 = [3835.673  3508.608  6931.577]

% Create a new matrix that includes locations where koala density exceeds the threshold
M2=[M(:,1:4),out2(:,1)];%[,1] = lon ; [,2]= lat ; [,3] = habitat_suitability ;	[,4] = mean_density ; [,5] = median value above the treshold


%% Decomposition function: calculate the proportion of total habitat suitability per location
% We use the proportion of habitat suitability to distribute the koala population across locations.
sumHBM = sum(M(:,3),"omitnan");decPop = (M(:,3).*100)/sumHBM;% Total habitat suitability (omit NaN values) + Percentage of total suitability per location
poptot=sum(M(:,4),"omitnan");decPop=[decPop,(decPop.*poptot)./100];% Total population of koalas (from the current data)


% Apply the future population stats from the demographic model (no intervention scenario, baseline)
futpoptot = stats(1);% Use baseline future population from stats(1)
decPopFut=[decPop(:,1),(decPop(:,1).*futpoptot)./100];% Redistribute future population based on habitat suitability
decPopFuttot = sum(decPopFut(:,2),"omitnan");% Ensure correct redistribution
PopFut = [M(:,1:3),decPopFut(:,2)];% Add the redistributed future population to coordinates and habitat suitability [,1] = lon ; [,2]= lat ; [,3] = habitat_suitability ;	[,4] = mean_density ;

% Check how many locations and how much of the koala population exceed the threshold of 70 per hectare in the future
outfut=PopFut(:,4)-tr;[r,c]=find(outfut<0);outfut2=outfut;% Find locations below threshold
for i=1:length(r);outfut2(r(i),c(i))=NaN;end;% Set future densities below threshold to NaN
statfut70=sum(outfut2,1,"omitnan"); % Calculate the total number of koalas above the threshold in the future

% Save data for plotting: [,1] = lon ; [,2] = lat ; [,3] = habitat suitability ; [,4] = mean density
% [,5] = median value above threshold (present day) ; [,6] = median value above threshold (future)
mapout=[M2,outfut2];save -ascii DensityMappingTreshold.txt mapout;

%% Run some statistics for proportion etc.
% Find the number of grid cells with koalas
id0=find(mapout(:,4)>0);lid0=length(id0);% Total number of grid cells with koalas
% Calculate the percentage of grid cells above the threshold for present day
TF1=isnan(mapout(:,5));id1=find(TF1==0);lid1=length(id1);pid1=(lid1.*100)/lid0;% Percentage of grid cells above threshold in present day
% Calculate the percentage of grid cells above the threshold for the future
TF2=isnan(mapout(:,6));id2=find(TF2==0);lid2=length(id2);pid2=(lid2.*100)/lid0;% Percentage of grid cells above threshold in the future


%%=========================================================================
%% STEP 2: IMPACT OF STERILISATION ON TOTAL POPULATION
%%=========================================================================
% This script provides a comprehensive visualisation of the impact of sterilisation on the koala population over time, comparing different intervention scenarios. Let me know if you need further clarification on any section!
% Showing the median future value in red as a function of sterilisation plan
% Showing in green where the population should be (future and present limits)

% Key Explanations:
%   Sterilisation Proportion (propster): The propster vector defines the sterilisation proportion from 0 to 0.9 in increments of 0.1.
%                                        The variable np holds the number of points in this vector.
%   Shading and Plotting Reference Lines:The plots show the impact of sterilisation on koala populations under different scenarios:
%                                        futpoptot: future total population
%                                        poptot: current population
%                                        futpoptot - statfut70: future population above the threshold (0.7 koalas per hectare)
%   Scenario 1 and Scenario 2:
%                                        Scenario 1: Only females are sterilised.
%                                        Scenario 2: Both females and their daughters are sterilised.
%                                        For each scenario, the median population values and confidence intervals (low and high CI) are plotted.
%                                        The vertical lines (xline) indicate the sterilisation proportion threshold for acceptable population control.
%   Color Schemes: brewermap is used to generate color palettes (clr for shading and clr1 for lines) for the plots.
%   Statistical Analysis: The script evaluates the impact of different sterilisation proportions on the koala population size, showing the median population size under different sterilisation scenarios.

% Proportion of individuals sterilised ranges from 0 to 0.9
propster=0:0.1:0.9;np=length(propster);% Number of points in the sterilisation proportion vector
clr=brewermap(9,'Pastel1');clr1=brewermap(12,'Paired');% Define color schemes using Brewermap for shading and lines
% Plot the limits of good and bad outcomes
figure(1),area(propster,repmat(futpoptot,np,1),'FaceColor',clr(1,:),'EdgeColor','none'),% Shading for future population total (good)
          hold on, area(propster,repmat(poptot,np,1),'FaceColor',clr(5,:),'EdgeColor','none'), % Shading for total population without intervention (bad)
          hold on, area(propster,repmat(futpoptot-statfut70,np,1),'FaceColor',clr(3,:),'EdgeColor','none'),% Shading for future population above threshold (bad)

%% Scenario 1: Female sterilisation only
% Change directory to where scenario 1 outputs are stored = to be modified as needed
cd '~/...'
% Load Scenario 1 median data and scale values by 2
f0med = readmatrix('Total_popFemale(Median)_scenario#1.csv');f0med(:,2:4)=f0med(:,2:4).*2;

% Define threshold for acceptable sterilisation proportion
x_value1 = 0.1205;
% Plot the sterilisation results for Scenario 1
figure(1),hold on, shade(propster, f0med(:,3), propster, f0med(:,4), 'FillType', [1, 2; 2, 1], 'FillColor', clr(2,:), 'FillAlpha', 0.6);% Shaded area for uncertainty (low to high CI)
          hold on, plot(propster, f0med(:,2),'Color',clr1(2,:),'LineWidth',1.5);% Plot median future population
          hold on, plot(propster, f0med(:,3),':k','Color',clr1(1,:),'LineWidth',1.5); % Plot lower CI
          hold on, plot(propster, f0med(:,4),':k','Color',clr1(1,:),'LineWidth',1.5); % Plot upper CI
          % Plot reference lines for future population, current population, and future population above threshold
          hold on, plot(propster,repmat(futpoptot,np,1),'Color',clr1(6,:),'LineWidth',1.5),% Future population total
          hold on, plot(propster,repmat(poptot,np,1),'Color',clr1(8,:),'LineWidth',1.5),% Total population without intervention
          hold on, plot(propster,repmat(futpoptot-statfut70,np,1),'Color',clr1(4,:),'LineWidth',1.5),% Future population above the threshold
          % Plot a vertical line indicating the sterilisation threshold (acceptable proportion)
          xline(x_value1, 'k--'); % Dashed black line at the sterilisation threshold
          axis square; xlab('proportion of sterilised individuals', 13);ylab('total population size (N)',13);% Configure plot appearance
          title('Koala Scenario #1 (median value) - treshold = 0.7','FontSize',15)

%% Scenario 2: Female and daughter sterilisation
% Change directory to where scenario 2 outputs are stored
figure(2),area(propster,repmat(futpoptot,np,1),'FaceColor',clr(1,:),'EdgeColor','none'),% Shading for future population total (good)
          hold on, area(propster,repmat(poptot,np,1),'FaceColor',clr(5,:),'EdgeColor','none'), % Shading for total population without intervention (bad)
          hold on, area(propster,repmat(futpoptot-statfut70,np,1),'FaceColor',clr(3,:),'EdgeColor','none'), % Shading for future population above threshold (bad)

cd '~/....'% Load Scenario 2 median data and scale values by 2
f1med = readmatrix('Total_popFemale(Median)_scenario#2.csv');f1med(:,2:4)=f1med(:,2:4).*2;
x_value2 = 0.084;% Define threshold for acceptable sterilisation proportion
% Plot the sterilisation results for Scenario 2
figure(2),hold on, shade(propster, f1med(:,3), propster, f1med(:,4), 'FillType', [1, 2; 2, 1], 'FillColor', clr(2,:), 'FillAlpha', 0.6);% Shaded area for uncertainty (low to high CI)
          hold on, plot(propster, f1med(:,2),'Color',clr1(2,:),'LineWidth',1.5);% Plot median future population
          hold on, plot(propster, f1med(:,3),':k','Color',clr1(1,:),'LineWidth',1.5); % Plot lower CI
          hold on, plot(propster, f1med(:,4),':k','Color',clr1(1,:),'LineWidth',1.5);% Plot upper CI
          % Plot reference lines for future population, current population, and future population above threshold
          hold on, plot(propster,repmat(futpoptot,np,1),'Color',clr1(6,:),'LineWidth',1.5),% Future population total
          hold on, plot(propster,repmat(poptot,np,1),'Color',clr1(8,:),'LineWidth',1.5),% Total population without intervention
          hold on, plot(propster,repmat(futpoptot-statfut70,np,1),'Color',clr1(4,:),'LineWidth',1.5),% Future population above the threshold
          % Plot a vertical line indicating the sterilisation threshold (acceptable proportion)
          xline(x_value2, 'k--'); % Dashed black line at the sterilisation threshold
          % Configure plot appearance
          axis square; xlab('proportion of sterilised individuals', 13);ylab('total population size (N)',13);
          title('Koala Scenario #2 (median value) - treshold = 0.7','FontSize',15)

%%=========================================================================
%% STEP 3: IDENTIFYING THE COSTS BASED ON POPULATION THRESHOLD
%%=========================================================================
% This script calculates and visualizes the costs of different sterilization strategies, providing a clear visual comparison between the two scenarios
% 1) Female only sterilization
% 2) Female + daughter sterilization
clear all, close all, clc,% Clear workspace, close all figures, and clear the command window

% Define the range for sterilization proportion (yg) and years (xg)
yg = 0:0.1:0.9; yg = yg';          % Sterilization proportions from 0 to 0.9 with steps of 0.1
xg = 2016:2040;                    % Years from 2016 to 2040
nyg = 0:0.01:0.9; yg = yg';        % Sterilization proportions from 0 to 0.9 with finer resolution (steps of 0.01)
nxg = 2016:0.1:2040;               % Years from 2016 to 2040 with finer resolution (steps of 0.1)
[X, Y] = meshgrid(nxg, nyg);       % Create a grid for interpolation over years and sterilization proportions

%% Cost for female-only sterilization (Scenario 1)
cd '~/...'% Change directory to the Scenario 1 results. To be modified as needed.
% Read in the sterilization cost data for Scenario 1
data2=readmatrix('MedianSterilisationCost_Years_(Median)_scenario#1.csv');data2=data2(:,2:end);[~,c2]=size(data2);% Skip the first column (if it's not part of the cost data) and get the number of columns in the cost data
mat2=[];for i=1:c2;mat2=[mat2;yg',repmat(xg(i),length(yg),1),data2(:,i)];end;% Reshape the cost data for interpolation = For each year, create a matrix of sterilization proportions and corresponding costs
rast2 = griddata(mat2(:,2),mat2(:,1),mat2(:,3),X,Y,'cubic');%Interpolate the data to fit the grid defined by X and Y.  Interpolating the cost data using cubic interpolation

y_value1 = 0.1205;  % Threshold sterilization proportion for Scenario 1


%% Cost for female + daughter sterilization (Scenario 2)
cd '~/...'% Change directory to the Scenario 1 results. To be modified as needed.
data1=readmatrix('MedianSterilisationCost_Years_(Median)_scenario#2.csv');data1=data1(:,2:end);[~,c]=size(data1);% Read in the sterilization cost data for Scenario 2, Skip the first column (if it's not part of the cost data) ansds Get the number of columns in the cost data
mat=[];for i=1:c;mat=[mat;yg',repmat(xg(i),length(yg),1),data1(:,i)];end;% Reshape the cost data for interpolation
rast = griddata(mat(:,2),mat(:,1),mat(:,3),X,Y,'cubic');% Interpolate the data to fit the grid defined by X and Y.  Interpolating the cost data using cubic interpolation

y_value2 = 0.084;% Threshold sterilization proportion for Scenario 2

%% Prepare the colorscale for plotting
mrast1 = max(max(rast));         % Maximum cost for Scenario 2
mrast2 = max(max(rast2));        % Maximum cost for Scenario 1
tmax = max([mrast1, mrast2]);    % Maximum value across both scenarios for setting color scale
clr = brewermap(13, 'YlGnBu');   % Colormap for visualization

%% Plot Scenario 2 (Female + Daughter sterilization)
figure(1),imagesc(flipud(rast)),% Plot the cost raster for Scenario 2 with the y-axis flipped
       hold on,[C,h]= contour(flipud(rast./1000),'LineColor','w','ShowText','on');  % Add contour lines for better visualization (divide cost by 1000 for readability)
       clabel(C,h,'FontSize',15,'Color','w'),ylim([0 90]);xlim([0 length(nxg)]);% Label the contour lines and set axes
       set(gca,'YTickLabel',flipud(yg'),'XTickLabel',[nxg(1),2020,2025,2030,2035]);
       tr1=title('Cost female + daughter');set(tr1,'Fontsize',14);clim([0,tmax]);colorbar;axis square;

%% Plot Scenario 1 (Female only sterilization)
figure(2),imagesc(flipud(rast2)), % Plot the cost raster for Scenario 1 with the y-axis flipped
       hold on,[C2,h2]= contour(flipud(rast2./1000),'LineColor','w','ShowText','on');% Add contour lines for better visualization (divide cost by 1000 for readability)
       clabel(C2,h2,'FontSize',15,'Color','w'),ylim([0 90]);xlim([0 length(nxg)]); % Label the contour lines and set axes
       set(gca,'YTickLabel',flipud(yg'),'XTickLabel',[nxg(1),2020,2025,2030,2035]);
       tr2=title('Cost female only');set(tr2,'Fontsize',14);clim([0,tmax]);colorbar;axis square;

%% Plot the thresholds for both scenarios       
figure(3),yline(y_value1, 'k--');          % Plot a horizontal line for the sterilization threshold in Scenario 1
       hold on;yline(y_value2, 'r--');          % Plot a horizontal line for the sterilization threshold in Scenario 2
       xlim([xg(1), xg(25)]);ylim([yg(1), yg(10)]);axis square;% set axes properties


%%=========================================================================
%% STEP 4: CALCULATION OF TOTAL COST OVER TIME
%%=========================================================================
% This section computes the total cost for each sterilisation scenario by finding the sterilisation proportion that matches the target values and summing the costs over time
% Explanation: 
%   vect1 and vect2: These vectors represent the absolute differences between the sterilisation proportions (Y(:,1)) and the target sterilisation proportions (y_value1 for Scenario 1 and y_value2 for Scenario 2). This is used to find the closest match in the grid.
%   id1 and id2: These are the indices of the rows in Y where the sterilisation proportion is closest to the desired value for each scenario. The cost calculation is based on these rows.
%   cost1 and cost2: These variables store the total cost over 25 years for each scenario. The cost values are summed across the corresponding sterilisation proportions, and any NaN values are ignored.

%% Scenario #1: Female-only sterilisation
vect1 = abs(Y(:,1) - y_value1); % Calculate the absolute difference between sterilisation proportions (Y column 1) and the target proportion (y_value1)
id1 = find(vect1 == min(vect1)); % Find the index (id1) of the row where the difference is minimized (i.e., closest to y_value1)

% Calculate the total cost over time for Scenario #1 by summing up the values in the cost matrix (rast2)
% corresponding to the identified sterilisation proportion (id1), ignoring NaN values
cost1 = sum(rast2(id1,:), 'omitnan'); % The total cost for female-only sterilisation over 25 years is 28,441,559.7 units

%% Scenario #2: Female + daughter sterilisation
vect2 = abs(Y(:,1) - y_value2); % Calculate the absolute difference between sterilisation proportions (Y column 1) and the target proportion (y_value2)
id2 = find(vect2 == min(vect2));% Find the index (id2) of the row where the difference is minimized (i.e., closest to y_value2)

% Calculate the total cost over time for Scenario #2 by summing up the values in the cost matrix (rast)
% corresponding to the identified sterilisation proportion (id2), ignoring NaN values
cost2 = sum(rast(id2,:), 'omitnan');% The total cost for female + daughter sterilisation over 25 years is 39,808,686.9 units








