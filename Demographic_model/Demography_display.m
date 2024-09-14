program Demography_display


%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% STEP #1: SCENARIO #1 = FEMALE ONLY
%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Key points:
%   Data Loading: The script loads median, low confidence interval (CI), and high CI values for the koala population in Scenario 1.
%   Scaling: The data is scaled to a new maximum value for plotting purposes.
%   Plotting: It generates two figures:
%       The first figure plots the median, high, and low values with shaded areas representing uncertainty.
%       The second figure repeats the process but uses the scaled data to fit a different range.
%       Color Schemes: Brewermap is used to define the color schemes for shading and the median lines.
%       Threshold Line: A vertical line at the threshold value is added to the plots to highlight the critical proportion of females sterilized.

%%PROPORTION ACTUAL POPULATION WITH STERILISATION
clear all, close all, clc,% Clear workspace, close all figures, and clear the command window
cd '~/....'% Change directory to the folder where unmanaged scenario data is stored. To be modified as needed

% Read the initial population data for Scenario 1 (median, low CI, and high CI values)
data1.md = readmatrix('PMin_initial_popFemale&Male(Median)_scenario1.csv');  % Median values for Scenario 1
data1.lw = readmatrix('PMin_initial_popFEmale&Male(LowCI)_scenario1.csv');  % Low confidence interval values for Scenario 1
data1.hg = readmatrix('PMin_initial_popFEmale&Male(HighCI)_scenario1.csv'); % High confidence interval values for Scenario 1


% Combine data into matrices for plotting
Med = [data1.lw(:,1), data1.lw(:,2), data1.md(:,2), data1.hg(:,2)];  % Matrix for median values (columns: time, low, median, high)
High = [data1.lw(:,1), data1.lw(:,4), data1.md(:,4), data1.hg(:,4)]; % Matrix for high values (columns: time, low, median, high)
Low = [data1.lw(:,1), data1.lw(:,3), data1.md(:,3), data1.hg(:,3)];  % Matrix for low values (columns: time, low, median, high)

% Define color schemes for plots
clr1 = brewermap(13, 'Greens');clr2 = brewermap(13, 'YlOrRd');

% Predefined maximum value from the median values
currentMax = max(Med(1,:));  % Get the maximum value from the first row of the median data
newMax = 1;  % Define the new maximum value for scaling
scalingFactor = newMax / currentMax;  % Calculate the scaling factor based on the current and new maximum

% Scale median, high, and low values to the new maximum
Med_scaled_to_newMax = Med * scalingFactor ;Med_scaled_to_newMax(:,1) = Med(:,1) ;% Scale the median values and  Retain the original time values in the first column
High_scaled_to_newMax = High * scalingFactor ;High_scaled_to_newMax(:,1) = High(:,1);% Scale the high values and  Retain the original time values in the first column
Low_scaled_to_newMax = Low * scalingFactor ;Low_scaled_to_newMax(:,1) = Low(:,1);% Scale the low values and  Retain the original time values in the first column

% Threshold for the proportion of females sterilised
x_value1 = 0.1205;  % This threshold represents the acceptable proportion of sterilised females

% Plot the data for Scenario 1 (Female Only)
figure,bg=area(High(:,1),[Low(:,2),High(:,4)-Low(:,2)]),set(bg(1),'Facecolor','none','EdgeColor','none'),set(bg(2),'Facecolor',clr1(2,:),'EdgeColor','none'),% Create shaded area between low and high values, Set the first section (low) with no color % Set the second section (high) with green color
       hold on,up=area(High(:,1),[High(:,2),High(:,4)-High(:,2)]),set(up(1),'Facecolor','none','EdgeColor','none'),set(up(2),'Facecolor',clr1(3,:),'EdgeColor','none'),% Plot additional areas for high values No color for the first section Green color for the second section
       hold on,lw=area(High(:,1),[Low(:,2),Low(:,4)-Low(:,2)]),set(lw(1),'Facecolor','none','EdgeColor','none'),set(lw(2),'Facecolor',clr1(3,:),'EdgeColor','none'),% Plot areas for low values No color for the first section Green color for the second section
       hold on,md=area(High(:,1),[Med(:,2),Med(:,4)-Med(:,2)]),set(md(1),'Facecolor','none','EdgeColor','none'),set(md(2),'Facecolor',clr2(4,:),'EdgeColor','none','FaceAlpha',0.3), % Plot the median values with shaded areas No color for the first section Set the median with orange color and transparency
       % Plot the median line
       hold on, plot(Med(:,1), Med(:,2), '-k', 'Color', clr2(6,:), 'Linewidth', 1.5);  % Plot the median line with red color
       % Add a vertical dashed line at the threshold value
       xline(x_value1, 'k--');  % 'k--' represents a black dashed line at the threshold
       % Set axis limits and labels
       ylim([0 1.3]),xlim([0 0.9]),axis square;ylab('proportion of founding population (male + female)',13),xlab('proportion female sterilised',13);
       title('Koala Scenario #1 (median value) - treshold = 0.7','FontSize',15)
       
% Plot the scaled values (Scenario 1 with scaling applied)       
figure,bg=area(High_scaled_to_newMax(:,1),[Low_scaled_to_newMax(:,2),High_scaled_to_newMax(:,4)-Low_scaled_to_newMax(:,2)]),set(bg(1),'Facecolor','none','EdgeColor','none'),set(bg(2),'Facecolor',clr1(2,:),'EdgeColor','none'),
       hold on,up=area(High(:,1),[High_scaled_to_newMax(:,2),High_scaled_to_newMax(:,4)-High_scaled_to_newMax(:,2)]),set(up(1),'Facecolor','none','EdgeColor','none'),set(up(2),'Facecolor',clr1(3,:),'EdgeColor','none'),
       hold on,lw=area(High(:,1),[Low_scaled_to_newMax(:,2),Low_scaled_to_newMax(:,4)-Low_scaled_to_newMax(:,2)]),set(lw(1),'Facecolor','none','EdgeColor','none'),set(lw(2),'Facecolor',clr1(3,:),'EdgeColor','none'),
       hold on,md=area(High(:,1),[Med_scaled_to_newMax(:,2),Med_scaled_to_newMax(:,4)-Med_scaled_to_newMax(:,2)]),set(md(1),'Facecolor','none','EdgeColor','none'),set(md(2),'Facecolor',clr2(4,:),'EdgeColor','none','FaceAlpha',0.3),
       %hold on,plot(Low(:,1),Low(:,2),'--g','Linewidth', 0.5);hold on,plot(High(:,1),High(:,2),'--g','Linewidth', 0.5);
       hold on,plot(Med_scaled_to_newMax(:,1),Med_scaled_to_newMax(:,2),'-k','Color', clr2(6,:),'Linewidth', 1.5);
       xline(x_value1, 'k--'); % 'k--' represents treshol to get into the ok box
       ylim([0 1.1]),xlim([0 0.9]),axis square;ylab('proportion of founding population (male + female)',13),xlab('proportion female sterilised',13);
       title('Koala Scenario #1 (median value) - treshold = 0.7','FontSize',15)
       


%%++++++++++++++++++++++++++++++++++++++++++++++++
%%PROPORTION ACTUAL POPULATION WITH STERILISATION
%Summary:
%   The script visualizes the proportion of sterilized females over time in Scenario 1, including median, low confidence interval (CI), and high CI values.
%   It uses shaded areas to represent the uncertainty between the low and high CI values, while a solid line represents the median values.
%   A vertical dashed line is drawn to indicate the sterilisation threshold at 0.1205.
%   The Brewermap color schemes are used to differentiate the CI and median values in the plot.

cd '~/....'% Change directory to the folder where unmanaged scenario data is stored. To be modified as needed

% Load data from CSV files for Scenario 1
% 'md' stores the median values, 'lw' stores the lower confidence interval (CI) values, and 'hg' stores the high CI values
data1.md = readmatrix('Proportion_Actual_Femalsteril(Median)_scenario1.csv');  % Load median values for Scenario 1
data1.lw = readmatrix('Proportion_Actual_Femalsteril(LowCI)_scenario1.csv');   % Load low CI values for Scenario 1
data1.hg = readmatrix('Proportion_Actual_Femalsteril(HighCI)_scenario1.csv');  % Load high CI values for Scenario 1

% Create matrices containing the time, low, median, and high values for plotting
Med = [data1.lw(:,1), data1.lw(:,2), data1.md(:,2), data1.hg(:,2)];  % Time, low, median, and high values for median data
High = [data1.lw(:,1), data1.lw(:,4), data1.md(:,4), data1.hg(:,4)]; % Time, low, median, and high values for high CI data
Low = [data1.lw(:,1), data1.lw(:,3), data1.md(:,3), data1.hg(:,3)];  % Time, low, median, and high values for low CI data

% Set color schemes for the plot
clr1=brewermap(15,'BrBg');clr2=brewermap(13,'YlOrRd'); % Brewermap color scheme: Brown-Blue-Green and Yellow-Orange-Red

% Define the threshold value for the proportion of sterilized females
x_value1 = 0.1205;  % This value represents the threshold for an acceptable proportion of sterilized females

% Plot the actual population with sterilisation using 'area' to create shaded areas
% Plot the high and low confidence intervals as shaded areas
figure,bg=area(High(:,1),[Low(:,2),High(:,4)-Low(:,2)]),set(bg(1),'Facecolor','none','EdgeColor','none'),set(bg(2),'Facecolor',clr1(7,:),'EdgeColor','none'), % Create shaded area between low and high values, No color for the first section (low CI), % Brown color for the second section (high CI)
       hold on,up=area(High(:,1),[High(:,2),High(:,4)-High(:,2)]),set(up(1),'Facecolor','none','EdgeColor','none'),set(up(2),'Facecolor',clr1(6,:),'EdgeColor','none','FaceAlpha',0.3),  % Create shaded area for high CI values, No color for the first section Light brown color with transparency for the second section
       hold on,lw=area(High(:,1),[Low(:,2),Low(:,4)-Low(:,2)]),set(lw(1),'Facecolor','none','EdgeColor','none'),set(lw(2),'Facecolor',clr1(6,:),'EdgeColor','none','FaceAlpha',0.3), % Create shaded area for low CI values No color for the first section Light brown color with transparency for the second section
       % Plot the median values as a solid line
       hold on, plot(Med(:,1), Med(:,2), '-k', 'Color', [128, 96, 77]./255, 'Linewidth', 1.5);  % Dark brown color for the median line
       % Add a vertical dashed line at the threshold value
       xline(x_value1, 'k--');  % Dashed black line at the sterilisation threshold (0.1205)
       % Set axis limits and labels
       ylim([0 350000]),xlim([0 0.9]),axis square;ylab('number of sterilised female',13),xlab('proportion female sterilised',13);
       title('Koala Scenario #1 (median value) - treshold = 0.7','FontSize',15)
  

%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% STEP #2: SCENARIO #2 = FEMALE + DAUGHTER ONLY
%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% This script handles the Scenario #2 for koala population analysis, focusing on female and daughter sterilisation. 
% It plots the proportion of the founding population sterilised and the number of sterilised females. Two different plots are created: 
% one showing scaled data for Scenario 2 and another showing the actual proportion of females sterilised.
% same comment as above apply here


%%PROPORTION ACTUAL POPULATION WITH STERILISATION,
cd '~/....'% Change directory to the folder where unmanaged scenario data is stored. To be modified as needed

% Read the initial population data for Scenario 2 (median, low CI, and high CI values)
data2.md = readmatrix('PMin_initial_popFemale&Male(Median)_scenario2.csv');  % Median values for Scenario 2
data2.lw = readmatrix('PMin_initial_popFEmale&Male(LowCI)_scenario2.csv');  % Low confidence interval values for Scenario 2
data2.hg = readmatrix('PMin_initial_popFEmale&Male(HighCI)_scenario2.csv'); % High confidence interval values for Scenario 2

% Combine data into matrices for plotting
Med2 = [data2.lw(:,1), data2.lw(:,2), data2.md(:,2), data2.hg(:,2)];  % Median data (columns: time, low, median, high)
High2 = [data2.lw(:,1), data2.lw(:,4), data2.md(:,4), data2.hg(:,4)]; % High data (columns: time, low, median, high)
Low2 = [data2.lw(:,1), data2.lw(:,3), data2.md(:,3), data2.hg(:,3)];  % Low data (columns: time, low, median, high)

% Predefined maximum value for scaling
currentMax2 = max(Med2(1,:));  % Get the maximum value from the first row of median data
scalingFactor2 = newMax / currentMax2;  % Calculate the scaling factor based on the new maximum

% Scale to the new maximum
Med_scaled_to_newMax2 = Med2 * scalingFactor2 ;Med_scaled_to_newMax2(:,1) = Med2(:,1) ;% Scale the median data and Retain the original time values
High_scaled_to_newMax2 = High2 * scalingFactor2 ;High_scaled_to_newMax2(:,1) = High2(:,1);% Scale the high data and Retain the original time values
Low_scaled_to_newMax2 = Low2 * scalingFactor2 ;Low_scaled_to_newMax2(:,1) = Low2(:,1);% Scale the low data and Retain the original time values

% Threshold for the proportion of females + daughters sterilised
x_value2 = 0.084;  % This threshold represents the acceptable proportion of females + daughters sterilised

% Plot the scaled data for Scenario 2 (Female + Daughter Only)
figure,bg2=area(High_scaled_to_newMax2(:,1),[Low_scaled_to_newMax2(:,2),High_scaled_to_newMax2(:,4)-Low_scaled_to_newMax2(:,2)]),set(bg2(1),'Facecolor','none','EdgeColor','none'),set(bg2(2),'Facecolor',clr1(2,:),'EdgeColor','none'),
       hold on,up2=area(High2(:,1),[High_scaled_to_newMax2(:,2),High_scaled_to_newMax2(:,4)-High_scaled_to_newMax2(:,2)]),set(up2(1),'Facecolor','none','EdgeColor','none'),set(up2(2),'Facecolor',clr1(3,:),'EdgeColor','none'),
       hold on,lw2=area(High2(:,1),[Low_scaled_to_newMax2(:,2),Low_scaled_to_newMax2(:,4)-Low_scaled_to_newMax2(:,2)]),set(lw2(1),'Facecolor','none','EdgeColor','none'),set(lw2(2),'Facecolor',clr1(3,:),'EdgeColor','none'),
       hold on,md2=area(High2(:,1),[Med_scaled_to_newMax2(:,2),Med_scaled_to_newMax2(:,4)-Med_scaled_to_newMax2(:,2)]),set(md2(1),'Facecolor','none','EdgeColor','none'),set(md2(2),'Facecolor',clr2(4,:),'EdgeColor','none','FaceAlpha',0.3),
       %hold on,plot(Low(:,1),Low(:,2),'--g','Linewidth', 0.5);hold on,plot(High(:,1),High(:,2),'--g','Linewidth', 0.5);
       hold on,plot(Med_scaled_to_newMax2(:,1),Med_scaled_to_newMax2(:,2),'-k','Color', clr2(6,:),'Linewidth', 1.5);
       xline(x_value2, 'k--'); % 'k--' represents treshol to get into the ok box
       ylim([0 1.1]),xlim([0 0.9]),axis square;ylab('proportion of founding population (male + female)',13),xlab('proportion female + daughter sterilised',13);
       title('Koala Scenario #2 (median value) - treshold = 0.7','FontSize',15)
       

%%++++++++++++++++++++++++++++++++++++++++++++++++
%%PROPORTION ACTUAL POPULATION WITH STERILISATION
cd '~/....'% Change directory to the folder where unmanaged scenario data is stored. To be modified as needed
data2.md=readmatrix('Proportion_Actual_Femalsteril(Median)_scenario2.csv');
data2.lw=readmatrix('Proportion_Actual_Femalsteril(LowCI)_scenario2.csv');
data2.hg=readmatrix('Proportion_Actual_Femalsteril(HighCI)_scenario2.csv');


Med2=[data2.lw(:,1),data2.lw(:,2),data2.md(:,2),data2.hg(:,2)];
High2=[data2.lw(:,1),data2.lw(:,4),data2.md(:,4),data2.hg(:,4)];
Low2=[data2.lw(:,1),data2.lw(:,3),data2.md(:,3),data2.hg(:,3)];
clr1=brewermap(15,'BrBg');bsv=-10;clr2=brewermap(13,'YlOrRd');

x_value2 = 0.084;

figure,bg2=area(High2(:,1),[Low2(:,2),High2(:,4)-Low2(:,2)]),set(bg2(1),'Facecolor','none','EdgeColor','none'),set(bg2(2),'Facecolor',clr1(7,:),'EdgeColor','none'),
       hold on,up2=area(High2(:,1),[High2(:,2),High2(:,4)-High2(:,2)]),set(up2(1),'Facecolor','none','EdgeColor','none'),set(up2(2),'Facecolor',clr1(6,:),'EdgeColor','none','FaceAlpha',0.3),
       hold on,lw2=area(High2(:,1),[Low2(:,2),Low2(:,4)-Low2(:,2)]),set(lw2(1),'Facecolor','none','EdgeColor','none'),set(lw2(2),'Facecolor',clr1(6,:),'EdgeColor','none','FaceAlpha',0.3),
       hold on,plot(Med2(:,1),Med2(:,2),'-k','Color', [128, 96, 77]./255,'Linewidth', 1.5);
       xline(x_value2, 'k--'); % 'k--' represents treshol to get into the ok box
       ylim([0 350000]),xlim([0 0.9]),axis square;ylab('number of sterilised female',13),xlab('proportion female sterilised',13);
       title('Koala Scenario #2 (median value) - treshold = 0.7','FontSize',15)
  



