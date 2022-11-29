%%% GPR in MATLAB // URI Phillip Parisi - Updated May 3, 2022
%%% Raw Decimation to create a downsampled dataset! Take every jth pt
clc, clear GPR_d hp metrics nnum X X_d Y Y_d, close all, format compact

%%% This code does a DECIMATION (downselection) to get SlimData

% Add gpr_functions to the path
% you can do this manually with addpath(.../filepath/gpr_functions) 
    % if below code does not work
dir_path = cd;
idcs = strfind(dir_path,'/');
func_dir = dir_path(1:idcs(end));
func_dir = strcat(func_dir,"gpr_functions");
addpath(func_dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD Random Training Dataset
%clc, clear all, format compact, close all

load("~/Data/approx_methods_datasets/matlab_gpr_datasets/GPR_dataset_Seamounts_20000pts_2D.mat")

%load("~/Data/approx_methods_datasets/matlab_gpr_datasets/A_random_training_data_gpr_40.mat")
%load("~/Data/approx_methods_datasets/matlab_gpr_datasets/B_random_training_data_gpr_5000.mat")
data_stats.nnum = length(X);
data_stats.beg = X(1);
data_stats.end = X(end);

% Create a Downsampled Dataset N -> M
b = 20; % keep every bth point
[X_d,Y_d] = decimationDownsample(X,Y,b);

clearvars -except nnum X Y X_d Y_d GPR data_stats

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SET Hyperparameters
hp.L = 50;
hp.sigma_p = 3;
hp.sigma_n = 2;
hp.kerneltype = 'exact';            % sparse or exact kernel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RUN GPR
disp('...runnning GPR...')

%%%%%% Execute GPR
%GPR = calcGPR(X,Y,hp,data_stats);
GPR_d = calcGPR(X_d,Y_d,hp,data_stats);
disp('...finished GPR!...')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Metrics + Visuals
if exist('GPR','var') && exist('GPR_d','var')
    metrics = calcComparisonMetrics(GPR,GPR_d);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Other



% figures of data obtained from this script are in 
% random_comparisons.m