% Kmeans clustering
% split data into m clusters, take centroids as datapoints

clc, clearvars, close all

% Raw Data
if ispc()
    addpath("..\..\data\")
else
    addpath("../../data/")
end

rawdata = readtable("wiggles_single_ping.csv");
data = rawdata(77:end-37,[2 3]);

% Params
dsample_percent = 0.1;

sample_vals = dsample_percent * height(data);

% Kmeans
kdata = table2array(data);
[idx,C] = kmeans(kdata,sample_vals);



% Plots
plot(data.Y,data.Z,'k.'), hold on
plot(C(:,1),C(:,2),'ro')
