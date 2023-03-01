% Uncertainty - constant vs. changing

%%% GPR in MATLAB // URI Phillip Parisi - Update June 2022
tic, clc, clearvars, close all, format compact

%%%% GUIDE TO USE
%%% .m files you need:
% this script, which is the main script
% gpr_functions folder (one directory above the mainscript)

% Add gpr_functions to the path (update path as needed!)
% you can do this manually with addpath(.../filepath/gpr_functions) 
    % or update this code based on your path
    % this adds the path one directory above mainscript
dir_path = cd;
idcs = strfind(dir_path,'\');
func_dir = dir_path(1:idcs(end));
func_dir = strcat(func_dir,"gpr_functions");
addpath(func_dir);

%%%% RUNNING & PARAMETERS TO TWEAK
% this script should be good to run out-of-the-box
% there will be randomly generated x-data (training data), plotted in blue
% the GPR data points are plotted in red w/ error bars (uncertainty)

% You can TUNE
% - Kernel Hyperparameters
% - nnum, number of points you want in the dataset


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SETUP

% Kernel Hyperparameters [not optimized/trained] & Noise
hp.L = 7;                   % lengthscale (high = smoother, low = noisier)
hp.sigma_p = 4;            % process noise (aka vertical scale, output scale)
hp.sigma_n = 1;            % sensor noise (used to create W)
hp.kerneltype = 'exact';     % 'exact' or 'sparse' approximate kernel

% Generate Training Data w/ Gaussian Noise (aka Raw Data)
X = (0:0.75:22)';       % vertical array, training X, uniform random
nnum = length(X);
noise.mu = 0; noise.sigma = hp.sigma_n;
%Y = 1*sin(2*pi/(0.5*nnum)*X) + (normrnd(noise.mu,noise.sigma,nnum,1));  % vertical array, training Y, sinusoidal + noise
Y = 10*sin(X/2) + normrnd(noise.mu,noise.sigma,nnum,1);

% Prediction Points (X_Star)
X_Star = [[(-10+X(1)):.1:(10+X(end))]'; X];            % vertical array, add training data for prediction X


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MATRIX CALCS

% Calculate V and Inv(V)                     % depends on training x-points only
W = (hp.sigma_n^2)*eye(nnum);                % Whitenoise (identity * sigmasquared)
W(4,4) = 10; W(7,7) = 10; W(10,10) = 10;
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
p2 = plot(sorted.X_Star,sorted.Y_Star_Hat,'k-','Linewidth',2); % prediction means
p3 = plot(X,Y,'ko','MarkerFaceColor','k','MarkerSize',4); % training data
xlabel('X Values'), ylabel('Y Values'), title('Gaussian Process Regression')
grid on, legend([p3 p2 p1],"Training Data","Prediction \mu","Prediction 2\sigma")

