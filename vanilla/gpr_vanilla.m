% Vanilla Execution of GPR by Rasmussen Intended

clc, clearvars, close all

% add path to GPR functoins
dir_path = cd;
idcs = strfind(dir_path,'\');
func_dir = dir_path(1:idcs(end));
func_dir = strcat(func_dir,"gpr_functions");
addpath(func_dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SETUP

%%% Hyperparameters
hp.L = 3.1;                   % lengthscale (high = smoother, low = noisier)
hp.sigma_p = 8.25;           % std dev -- process noise (aka vertical scale, output scale)
hp.sigma_n = 2;            % std dev -- sensor noise (used to create W), not really a hp but kept here as a param
hp.kerneltype = 'exact';     % 'exact' or 'sparse' approximate kernel


%%% Generate Training Data w/ Gaussian Noise (aka Raw Data)

% training data noise
noise.mu = 0; noise.sigma = hp.sigma_n; 

% training vectors (vertical arrays)
X = (0:0.5:10)';       
Y_raw = 10*sin(X/2);
Y = Y_raw + normrnd(noise.mu,noise.sigma,length(X),1);  % add noise

% prediction points (X_Star)
X_Star = (-10 + X(1)) : 0.1 : (10 + X(end)); % new points
X_Star = [X_Star'; X]; % also predict at training points
X_Star = sortrows(X_Star,'ascend');

hp.nnum = length(X); % store # of vals


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MATRIX CALCS

% Calculate V and Inv(V)                     % depends on training x-points only
W = (hp.sigma_n^2)*eye(hp.nnum);             % Whitenoise (identity * sigmasquared)
V = K_Function(X,X,hp) + W;                  % Calculate Covariance Matrix using Kernel

% Generate K Parameters
K_Star = K_Function(X_Star,X,hp);            % Calculate K_Star for New Point(s)
K_StarStar = K_Function(X_Star,X_Star,hp);   % Calculate K_StarStar for New Point(s)

% Calculate Predictions!                                    % Finally bring in the training y-points here
Y_Star_Mean = K_Star * inv(V) * Y;                  % Mean Predictions (mean of Gaussians)
CapSigma_Star = K_StarStar - K_Star * inv(V) * K_Star' ; % Variance Predictions (prediction covariance matrix)
Y_Star_Var = diag(CapSigma_Star);                           % The diagonals store the variances we want!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOG MARGINAL LIKELIHOOD

LML = -1/2*Y'*inv(V)*Y - 1/2*log(det(V)) - 1/2*length(X)*log(2*pi)
disp('one iteration completed')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS & OUTPUTS

%%% Plots

% raw data points (with noise in Y)
figure(1)
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

%fill([X_Star fliplr(X_Star)], [lower_conf_bound fliplr(upper_conf_bound)], 'c')
%fill(X_Star,upper_conf_bound,'c','facealpha',0.3)
%{
%%% Plot GPR
% Organize Data to Plot It (created new obj
sortobj = [X_Star, Y_Star_Hat, Y_Star_Var];
sortobj = sortrows(sortobj);
sorted.X_Star = sortobj(:,1); sorted.Y_Star_Hat = sortobj(:,2); sorted.Y_Star_Var = sortobj(:,3);

% Bounded Plot (Training Data + Predictions + 2sigma Upper and Lower Bound
figure
p1 = plot(sorted.X_Star,sorted.Y_Star_Hat + 2*sqrt(sorted.Y_Star_Var),'r','LineWidth',2); hold on %upper bound
plot(sorted.X_Star,sorted.Y_Star_Hat - 2*sqrt(sorted.Y_Star_Var),'r','LineWidth',2); % lower bound
p2 = plot(sorted.X_Star,sorted.Y_Star_Hat,'k-','Linewidth',2); % prediction means
p3 = plot(X,Y,'ko','MarkerFaceColor','k','MarkerSize',4); % training data
xlabel('X Values'), ylabel('Y Values'), title('Gaussian Process Regression')
grid on, legend([p3 p2 p1],"Training Data","Prediction \mu","Prediction 2\sigma")
%}
