%%% GPR in MATLAB // URI Phillip Parisi - Updated May 3, 2022

%%% How to use this code
% 1. add function path
% 2. load data
% 2. run standard GPR (no downsampling)

clc, clear all, close all, format compact

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Add File Paths
dir_path = cd;
idcs = strfind(dir_path,'/');
main_dir = dir_path(1:idcs(end));
func_dir = [strcat(main_dir,"gpr_functions"), strcat(main_dir,"other_functions")];
addpath(func_dir(1)), addpath(func_dir(2))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD Random Training Dataset

% Load Single Ping Data
ping_filename = "wiggles_single_ping.csv";
data = readtable(strcat(main_dir,"/data/",ping_filename));
training.x = table2array(data(:,1));
training.y = table2array(data(:,2));
training.z = table2array(data(:,3));

training.x = training.x - mean(training.x); % de-mean x
training.y = training.y - mean(training.y); % de-mean y
training.z = training.z - mean(training.z); % de-mean z
training.npts = length(training.x);


disp(strcat("...training data with ", num2str(training.npts), " datapoints loaded..."))

% if you want new random data, scroll to bottom of this script

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RUN GPR on the Data
close all
disp('...runnning GPR...'), tic

X = training.x; Y = training.z; nnum = training.npts;
X_beg = X(1); X_end = X(end);
%%%%%%% Matric Calculations

% Kernel Hyperparameters [not optimized/trained] & Noise
hp.L = 0.3;                  % lengthscale (high = smoother, low = noisier)
hp.sigma_p = 0.23;            % process noise (aka vertical scale, output scale)
hp.sigma_n = 0.04;            % sensor noise (used to create W)
hp.kerneltype = 'exact';     % 'exact' or 'sparse' approximate kernel


% Prediction Points (X_Star)
X_Star = [[(-1+X_beg):.02:(1+X_end)]'; X];            % vertical array, add training data for prediction X


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MATRIX CALCS

% Calculate V and Inv(V)                     % depends on training x-points only
W = (hp.sigma_n^2)*eye(nnum);                % Whitenoise (identity * sigmasquared)
V = K_Function(X,X,hp) + W;                  % Calculate Covariance Matrix using Kernel

% Generate K Parameters
K_Star = K_Function(X_Star,X,hp);            % Calculate K_Star for New Point(s)
K_StarStar = K_Function(X_Star,X_Star,hp);   % Calculate K_StarStar for New Point(s)

% Cholesky Decomposition
L = chol(V,'lower');                         % Lower triangular cholesky factor

% Calculate Predictions!                                    % Finally bring in the training y-points here
Y_Star_Hat = K_Star * CholeskySolve(L,Y);                   % Mean Predictions (mean of Gaussians)
CapSigma_Star = K_StarStar-K_Star*CholeskySolve(L,K_Star'); % Variance Predictions (prediction covariance matrix)
Y_Star_Var = diag(CapSigma_Star);                           % The diagonals store the variances we want!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOG MARGINAL LIKELIHOOD

% How good is our fit? Use this to tune hyperparameters
LML = calcLML(L,Y,nnum);
AlgoTime = toc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS & OUTPUTS

% Output LML
fprintf('AlgoTime = %1.2f.\n',AlgoTime)
fprintf('Log Marginal Likelihood is %1.1f. Tune hyperparams for better fit.\n',LML)

%%% Plot GPR
% Organize Data to Plot It (created new obj
sortobj = [X_Star, Y_Star_Hat, Y_Star_Var];
sortobj = sortrows(sortobj);
sorted.X_Star = sortobj(:,1); sorted.Y_Star_Hat = sortobj(:,2); sorted.Y_Star_Var = sortobj(:,3);

% Bounded Plot (Training Data + Predictions + 2sigma Upper and Lower Bound
figure
p1 = plot(sorted.X_Star,sorted.Y_Star_Hat + 2*sqrt(sorted.Y_Star_Var),'r','LineWidth',2); hold on %upper bound
plot(sorted.X_Star,sorted.Y_Star_Hat - 2*sqrt(sorted.Y_Star_Var),'r','LineWidth',2); % lower bound
p2 = plot(sorted.X_Star,sorted.Y_Star_Hat,'r--','Linewidth',2); % prediction means
p3 = plot(X,Y,'bo','MarkerFaceColor','b','MarkerSize',3); % training data
xlabel('Position on Seafloor'), ylabel('Depth'), title('Gaussian Process Regression')
grid on, legend([p3 p2 p1],"Training Data","Prediction \mu","Prediction 2\sigma")

axis equal













%%
%%%%%% Calculate Metrics

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