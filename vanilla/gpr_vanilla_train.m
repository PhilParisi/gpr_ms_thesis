% Vanilla Execution of GPR by Rasmussen Intended
% meant to be a walk through of running GPR
% run each section by section and understand what is happening

clc, clearvars, close all

% add path to GPR functions
dir_path = cd;
idcs = strfind(dir_path,'\');
func_dir = dir_path(1:idcs(end));
func_dir = strcat(func_dir,"gpr_functions");
addpath(func_dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SETUP

% Hyperparams to loop over in optimization
L = 0.1:0.5:20;
sigma_p = 0.01:0.2:3;

%%% Generate Training Data w/ Gaussian Noise (aka Raw Data)

% training data noise
noise.mu = 0; noise.sigma = 2; 

% training data vectors (vertical arrays)
X = (0:0.5:10)';       
Y_raw = 10*sin(X/2);
Y = Y_raw + normrnd(noise.mu,noise.sigma,length(X),1);  % add noise

% prediction points (X_Star)
X_Star = (-5 + X(1)) : 0.1 : (5 + X(end)); % new points
X_Star = [X_Star'; X]; % also predict at training points
X_Star = sortrows(X_Star,'ascend');

% start a hyperparameters structure
hp.nnum = length(X); % store # of vals



%%%% OPTIMIZATION
disp('running optimization routine, I will tell you when Im done')

for i = 1:length(L)
    for j = 1:length(sigma_p)

% set hyperparams
hp.L = L(i);
hp.sigma_p = sigma_p(j);
hp.sigma_n = noise.sigma;
hp.kerneltype = "exact";


% GPR equations

% calculate V and inv(V)                     % depends on training x-points only
W = (hp.sigma_n^2)*eye(hp.nnum);             % Whitenoise (identity * sigmasquared)
V = K_Function(X,X,hp) + W;                  % Calculate Covariance Matrix using Kernel

% generate K matrices
K_Star = K_Function(X_Star,X,hp);            % Calculate K_Star for New Point(s)
K_StarStar = K_Function(X_Star,X_Star,hp);   % Calculate K_StarStar for New Point(s)

% calculate predictions!                                    % Finally bring in the training y-points here
Y_Star_Mean = K_Star * inv(V) * Y;                  % Mean Predictions (mean of Gaussians)
CapSigma_Star = K_StarStar - K_Star * inv(V) * K_Star' ; % Variance Predictions (prediction covariance matrix)
Y_Star_Var = diag(CapSigma_Star);                           % The diagonals store the variances we want!


% calculate log marginal likelihood
LML = -1/2*Y'*inv(V)*Y - 1/2*log(det(V)) - 1/2*length(X)*log(2*pi);
LML_map(i,j) = LML;

    end
end

disp('done calculating hyperparam combos, run next section!')


%% View Optimization Results
close all

% create LML surface plot
figure()
[L_grid,sigma_p_grid] = meshgrid(L,sigma_p);
surf(L_grid',sigma_p_grid',LML_map)
ylabel('process noise'),xlabel('lengthscale'),zlabel('LML')

% create LML contour plot
figure()
contour(L', sigma_p', LML_map')
xlabel('lengthscale'),ylabel('process noise')
title('LML Likelihood Contour'), grid on

% find maximum value and corresponding hyperparams
[max_value, max_index] = max(LML_map(:));
[max_row, max_col] = find(LML_map == max_value);

% display to user the best hyperparams
fprintf('Best lengthscale: %d\n', L(max_row));
fprintf('Best sigma_p: %d\n', sigma_p(max_col));

% other output messages
disp('now run the next section using these hyperparams')


%% Trained GPR
% yes, we're running the same GPR equations again


% set hyperparams
hp.L = L(max_row);
hp.sigma_p = sigma_p(max_col);
hp.sigma_n = noise.sigma;
hp.kerneltype = "exact";

% GPR equations

% calculate V and inv(V)                     % depends on training x-points only
W = (hp.sigma_n^2)*eye(hp.nnum);             % Whitenoise (identity * sigmasquared)
V = K_Function(X,X,hp) + W;                  % Calculate Covariance Matrix using Kernel

% generate K matrices
K_Star = K_Function(X_Star,X,hp);            % Calculate K_Star for New Point(s)
K_StarStar = K_Function(X_Star,X_Star,hp);   % Calculate K_StarStar for New Point(s)

% calculate predictions!                                    % Finally bring in the training y-points here
Y_Star_Mean = K_Star * inv(V) * Y;                  % Mean Predictions (mean of Gaussians)
CapSigma_Star = K_StarStar - K_Star * inv(V) * K_Star' ; % Variance Predictions (prediction covariance matrix)
Y_Star_Var = diag(CapSigma_Star);                           % The diagonals store the variances we want!


% calculate log marginal likelihood
LML = -1/2*Y'*inv(V)*Y - 1/2*log(det(V)) - 1/2*length(X)*log(2*pi);


%%% Plots

% raw data points (with noise in Y)
figure()
plot(X,Y,'bo','markerfacecolor','b'), hold on

% raw data points / true curve (with no noise, the true function)
plot(X,Y_raw,'k-')

% raw data error bars (2*std = 95% confidence)
input_std_devs = sqrt(diag(W)); % W is 
input_error_bars = 2*input_std_devs; % 2*std dev is 95% confidence 
errorbar(X,Y,input_error_bars,'ob')

% prediction means
plot(X_Star,Y_Star_Mean, 'r')

% prediction confidence bounds (2*std = 95% confidence)
upper_conf_bound = Y_Star_Mean + 2*sqrt(Y_Star_Var);
lower_conf_bound = Y_Star_Mean - 2*sqrt(Y_Star_Var);
plot(X_Star,upper_conf_bound,'r')
plot(X_Star,lower_conf_bound,'r')