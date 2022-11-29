%%% GPR in MATLAB // URI Phillip Parisi - Updated May 3, 2022
%%% Approach from Dr. Kristopher Krasnosky

%%% Version Control Notes
% Working with approximate methods now
% Added in Cholesky option to calcGPR()

%%% How to use this code
% 1. either load a saved dataset OR calculate new dataset
% 2. run GPR section


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD Random Training Dataset
clc, clear all, close all, format compact
addpath("gpr_functions/")

%load("~/Data/approx_methods_datasets/matlab_gpr_datasets/A_random_training_data_gpr_40.mat")
load("~/Data/approx_methods_datasets/matlab_gpr_datasets/B_random_training_data_gpr_5000.mat")
nnum = length(X);

% Create a Downsampled Dataset N -> M
b = 10; % keep every jth point
[X_d,Y_d] = rawDownsample(X,Y,b);

disp(strcat("...training data with ", num2str(nnum), " datapoints loaded..."))
disp(strcat("...downsampled: kept every 1 in ", num2str(b), " points..."))
clearvars -except nnum X Y X_d Y_d

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NEW Random Training Dataset
clc, clear all, close all, format compact

% Training Data + Noise
nnum = 5000; 
X_beg = -nnum; X_end = nnum;                        
X = (X_end - X_beg)*rand(nnum,1) + X_beg;           % vertical array, training X, uniform random
X = sort(X);                                        % put in ascending order

Y = 3*sin(2*pi/40*X) + (rand(nnum,1)*1 - 0.5);      % vertical array, training Y, sinusoidal + noise

% Create a Downsampled Dataset N -> M
b = 2; % keep every jth point
[X_d,Y_d] = rawDownsample(X,Y,b);

disp(strcat("...new training data generated with ",num2str(nnum)," points..."))
disp(strcat("...downsampled: kept every 1 in ", num2str(b), " points..."))
clearvars -except nnum X Y X_d Y_d

%% RUN GPR
close all
disp('...runnning GPR...')

%%%%%% Execute GPR
% Use either cholesky decomp -->    calcGPR_chol()
% or regular inverse with inv() --> calcGPR()
GPR = calcGPR(X,Y,nnum,'chol');
GPR_d = calcGPR(X_d,Y_d,nnum,'chol');
disp('...finished GPR!...')

%%%%%% Calculate Metrics
% each method against TRAINING DATA? GRIDDED PRODUCT?

%%% Sum of Squared Differences - between exact and approx
Metrics.SSD = sum((GPR.Y_Star_Hat - GPR_d.Y_Star_Hat).^2);

%%% R-Squared - between exact and approx
Metrics.SST = sum((GPR.Y_Star_Hat - ones(GPR.npredict,1)*mean(GPR.Y_Star_Hat)).^2);
Metrics.Rsqr = 1 - Metrics.SSD / Metrics.SST;

%%% Calculate Speed Comparison
Metrics.PerFaster = GPR.compute_time / GPR_d.compute_time;

%%% Output Results to Command Window
fprintf("SSD (Exact-Approx): %0.3f\n",Metrics.SSD)
fprintf("RSqr (Exact-Approx): %0.3f\n",Metrics.Rsqr)
fprintf("Times Faster (Exact-Approx): %0.3f\n",Metrics.PerFaster)

%%% Calculate Speed Comparison
Metrics.PerFaster = GPR.compute_time / GPR_d.compute_time;

%%%%%% Data Visualization

% Combined Plot
figure
subplot(2,1,1)
plot(GPR.X_Star,GPR.Y_Star_Hat,'ro','MarkerFaceColor','r','MarkerEdgeColor','k'), hold on
errorbar(GPR.X_Star,GPR.Y_Star_Hat,GPR.Sigma_Star,'r.','LineWidth',2)
plot(GPR.X,GPR.Y,'bo','MarkerFaceColor','b','MarkerSize',8), grid on, hold on
xlabel('X Values'), ylabel('Y Values')
legend('Predicted', 'Uncertainty','Raw Data')
title(GPR.title_str)

subplot(2,1,2)
plot(GPR_d.X_Star,GPR_d.Y_Star_Hat,'ro','MarkerFaceColor','r','MarkerEdgeColor','k'), hold on
errorbar(GPR_d.X_Star,GPR_d.Y_Star_Hat,GPR_d.Sigma_Star,'r.','LineWidth',2)
plot(GPR_d.X,GPR_d.Y,'bo','MarkerFaceColor','b','MarkerSize',8), grid on, hold on
xlabel('X Values'), ylabel('Y Values')
legend('Predicted', 'Uncertainty','Raw Data')
title(GPR_d.title_str)

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