% Hyperpameters Visualizer
% the csv_gpr_gui outputs a csv with all the hp's tested
% this script opens the csv and find the max LML

clc, clearvars, close all

if ispc()
    addpath('..\results\LML')
else
    addpath('../results/LML/')
end


% load the csv
data = readtable("exact_100_lml.csv");
%data = readtable("exact_100_lml.csv");
% remove rows that have NAN LML
%to_remove = isnan(data.lml)';
%data = data(~to_remove,:);

% output the max lml value
[lml_max, ind_max] = max(data.lml);
fprintf("max lml is %d at length=%d and process=%d \n",lml_max,data.length_scale(ind_max),data.process_noise(ind_max));


%%
% make a scatter plot
close all
scatter3(data.length_scale,data.process_noise,data.lml)
xlabel('lengthscale'),ylabel('process_noise'),zlabel('lml')