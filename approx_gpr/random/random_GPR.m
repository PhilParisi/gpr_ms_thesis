%%% GPR in MATLAB // URI Phillip Parisi - Updated June 7, 2022
clc, clear all, close all
%% UNIFORM RANDOM DOWNSAMPLE RAW DATA INTO SLIM DATA
clc, clear GPR_d hp metrics nnum X X_d Y Y_d, close all, format compact

%%% This code does a RANDOM downselection to get SlimData
% Given a # of desired pts in the SlimData, indices are chosen by a uniform
%   random variable
% A random subset should have similar properties to the larger dataset

% Add gpr_functions to the path
% you can do this manually with addpath(.../filepath/gpr_functions) 
    % if below code does not work
dir_path = cd;
idcs = strfind(dir_path,'/');
func_dir = dir_path(1:idcs(end));
func_dir = strcat(func_dir,"gpr_functions");
addpath(func_dir);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD Training Dataset
clc, format compact, close all

load("~/Data/approx_methods_datasets/matlab_gpr_datasets/GPR_dataset_Seamounts_20000pts_2D.mat")
%load("~/Data/approx_methods_datasets/matlab_gpr_datasets/B_random_training_data_gpr_5000.mat")

data_stats.nnum = length(X);
data_stats.beg = X(1);
data_stats.end = X(end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RANDOM Downsample
slim_pts = round(0.05*data_stats.nnum); % #pts in slim data
[X_d,Y_d] = randomDownsample(X,Y,slim_pts);
%clearvars -except nnum X Y X_d Y_d


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SET Hyperparameters
hp.L = 1000;
hp.sigma_p = 3;
hp.sigma_n = 0.6;
hp.kerneltype = 'sparse';        %exact or sparse kernel


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






%% Data from Speed Trials [Figure Code]
%clc, close all 

% TRIALS RUN ON MAY 3rd 2022
% MATLAB Without Cholesky Decomposition
dpts =  [2500,  1667,   1250,   1000,   834,    715,    625,    556,    500];
R2 =    [.997,  .943,   .886,   .814,   .741,   .667,   .586,   .512,   .503];
fast =  [2.479, 3.672,  3.967,  4.684,  4.866,  5.362,  5.673,  5.908,  6.059]; % updated May 3
% With Cholesky Decomposition
fast_c= [2.772, 4.356,  5.037,  5.675,  6.570,  6.554,  6.690,  7.121,  7.313]; % updated May 3


figure
yyaxis left
plot(dpts,fast,'--ob','MarkerFaceColor','b'), hold on
plot(dpts,fast_c,'--sm','MarkerFaceColor','m')
ylabel('Times Faster')
ylim([2, 8])
hold on
yyaxis right
plot(dpts,R2,'--vr','MarkerFaceColor','r')
ylabel('R-Squared Fit')
grid on
set(gca,'xdir','rev')
title('Comparing Different Raw Downsampled points, Exact had 5000pts')
xlabel('Downsampled dpts')
legend('Inverse','w/ Cholesky')