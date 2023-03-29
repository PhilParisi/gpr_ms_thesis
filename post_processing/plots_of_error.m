% The Script that Makes all the Error Plots
% the data below is all from training_hp.xlsx

clc, clearvars, close all

% EXACT
exact.time = 92.06;

% AVERAGE DOWNSAMPLE
average.dsample = [50, 33, 25, 20, 10];
average.RMSE = [0.0580, 0.0683, 0.0847, 0.1047, 0.1663];
average.MAE = [0.0355, 0.0412, 0.0460, 0.0531, 0.0953];
average.rawtime = [44.53, 28.54, 22.69, 18.71, 10.74]; 
average.tiletime = [0.128, 0.093, 0.075, 0.063, 0.044];

for i = 1:length(average.dsample)
    plot(average.rawtime(i), average.RMSE(i),'bo','markerfacecolor','b')
    hold on
end

xlabel('time (secs)'),ylabel('RMSE'), grid on
hold off
title('Point Averaging Downsample Results')

for i = 1:length(average.dsample)
    gtext(char(strcat(num2str(average.dsample(i)),'%')))
end
gtext('Exact took 92 seconds')



%%
clc, clearvars, close all 

%%%%%% HYBRID DOWNSAMPLE
hybrid.dsample = [90, 80, 70, 60, 50, 40, 30, 20, 10];
hybrid.RMSE = [0.0851, 0.1047, 0.0727, 0.0851, 0.0979, 0.0964, 0.0864, 0.1162, 0.1659];
hybrid.MAE = [0.0522, 0.0595, 0.0436, 0.0509, 0.0576, 0.0586, 0.0526, 0.0676, 0.1009];
hybrid.rawtime = [110.15, 89.53, 73.28, 65.94, 54.77, 43.44, 30.46, 20.50, 10.83]; 
hybrid.tiletime = [0.372, 0.309, 0.245, 0.227, 0.202, 0.161, 0.110, 0.083, 0.056];

plot(hybrid.tiletime, hybrid.RMSE,'bo','markerfacecolor','b')


%%

