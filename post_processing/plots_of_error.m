% The Script that Makes all the Error Plots
% the data below is all from training_hp.xlsx

clc, clearvars, close all

% EXACT
exact.time = 92.06;

% AVERAGE DOWNSAMPLE
average.dsample = [50, 33, 25, 20, 10];
average.RMSE = [0.0580, 0.0683, 0.0847, 0.1047, 0.1663];
average.MAE = [0.0355, 0.0412, 0.0460, 0.0531, 0.0953];
average.time = [44.53, 28.54, 22.69, 18.71, 10.74]; 

for i = 1:length(average.dsample)
    plot(average.time(i), average.RMSE(i),'bo','markerfacecolor','b')
    hold on
end

xlabel('time (secs)'),ylabel('RMSE'), grid on
hold off
title('Point Averaging Downsample Results')

for i = 1:length(average.dsample)
    gtext(char(strcat(num2str(average.dsample(i)),'%')))
end
gtext('Exact took 92 seconds')