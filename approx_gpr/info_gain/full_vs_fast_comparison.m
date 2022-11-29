% Comparing Fast GPR (SlimData) vs. Full GPR  (RawData)
% Programmer: Parisi
% Date: 31May2022

% GOAL
    % This script aims to generate figures of relative comparison between
    % different GPR methods (full vs fast)

%%
addpath("gpr_functions/")

%%
%%%%%%%%% Bring in Raw Data
clc, clear all, format compact, %close all

%seamounts
load("~/Data/approx_methods_datasets/matlab_gpr_datasets/GPR_dataset_Seamounts_20001pts_2D.mat")

%load("~/Data/approx_methods_datasets/matlab_gpr_datasets/A_random_training_data_gpr_40.mat")
%load("~/Data/approx_methods_datasets/matlab_gpr_datasets/B_random_training_data_gpr_5000.mat")
nnum = length(X);

disp(strcat("...training data with ", num2str(nnum), " datapoints loaded..."))
clearvars -except nnum X Y

% X and Y must be VERTICAL!

%% SELECT MULTIPLE POINTS (GROUP OF INCLUSION POINTS)

clc, clearvars -except nnum X Y, close all
%Downsampling as Described in Paper

% total points wanted in slim data [threshold for ending algo]
totalpts_slimdata = round(0.25*nnum);

%%%%%%%%%  Hyperparamters
sigma = 0.8;
sigmasq = sigma^2;
% lengthscale and process noise currently stored in kernel function

%%%%%%%%%  Indexing Variables
I = []; % Active Set of Indices to Include (starts empty), w/ d values
    % add to this list with --> I(end+1) = value
R = 1:nnum; % Remaining Indices not in Active Set, I
    % remove from this list with --> R(index) = []

%%%%%%%%%  Initialize a Random Value into I to start
%inclusion_ind = 27;             %fixed
inclusion_ind = randi(nnum,1);  %random

I = [I inclusion_ind];
R(R==inclusion_ind) = [];


%%%%%%%%%  CURRENT TERMS

% beginning terms
K_I_I = K_Function(X(I),X(I));      % Self Covariance, I vs. I, dxd
K_I_dot = K_Function(X(I),X(:));    % Joint Covariance, I vs. All, dxn
P_I = inv(K_I_I) * K_I_dot;         % !?!?INVERSE!?!? dxd * dxn

% cholesky decomp of K_I_I
L = chol(K_I_I,'lower');        % dxd    % Lower Triangular
V = inv(L)*K_I_dot;             % dxn
%Ident_I = eye(length(I));       % dxd    % Identity Matrix
M = sigmasq*eye(length(I)) + V*V';     % dxd    % Noisy Covariance maybe? 

% cholesky decomp of M
L_M = chol(M,'lower');          % dxd
Beta = inv(L_M)*V*Y;            % dx1

% p and q
p = diag(V'*V);                 % nx1
q = diag(V'*inv(M)*V);          % nx1

% mu
mu = V'*inv(L_M')*Beta;         % ANOTHER INVERSE???


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ITERATION MODEL DEVELOPMENT


%%%%%%%%%%  OUTER LOOP
while length(I) < totalpts_slimdata         % max pts in Slim Data



    %%%%%%%%%%  INFORMATION GAIN for Each Inclusion Candidate
    % go thru each i in R, and calculate info gain
    info_gain = [];

    for R_ind = R       % go thru all R indices

        % APPROXIMATE INFORMATION GAIN
        p_i = p(R_ind);
        q_i = q(R_ind);

        % This is becoming imaginary when p_i > sigmasq
        l_i = sqrt(round(K_Function(X(R_ind),X(R_ind)) - p_i, 5)); 

        epsilon_i = 1/((sigma/l_i)^2 + 1 - q_i);
        kappa_i = epsilon_i*(1+2*(sigma/l_i)^2);

        % info gain equation
        info_gain(R_ind,1) = -log(sigma/l_i) + ... % NATURAL LOG?
            -0.5*(log(epsilon_i) + ...
            epsilon_i*(1-kappa_i)*(sigma^-2)*((Y(R_ind)-mu(R_ind))^2) + ...
            -kappa_i + 2);
            % recall log is natural log, log10 is log

    end % output is a info_gain column vector


    %%%%%%%%%%  INCLUSION into MODEL
    % Inclusion Point
    [inclusion_val, inclusion_ind] = max(info_gain);


    %%%%% UPDATE SETS I, R
    I = [I inclusion_ind];      % add inclusion candidate to I
    R(R==inclusion_ind) = [];   % remove inclusions candidate from R


    %%%%% Update Terms to include the newly selected candidate, e.g. I --> I_prime

    % UPDATE L (dxd) to L_prime (d+1 x d+1), recall L is lower triangular
    L(end+1,:) = (V(:,inclusion_ind))'; % new row
    %this term can sometimes go negative!!!
    L(end,end+1) = sqrt( K_Function(X(inclusion_ind),X(inclusion_ind)) - p(inclusion_ind) ); % bottom rt scalar
    
    li = L(end,end); % bottom rt corner of Lprime
    vi = V(:,inclusion_ind); %dx1 col vector
    

    % UPDATE V (dxn) to V_prime (d+1 x n, add new row)
    v = 1/li * (K_Function(X,X(inclusion_ind)) - V'*vi); %nx1 col vector
    V_prime = [V; v'];


    p = p + v.^2; % is this the right approach?

    % UPDATE L_M (dxd to d+1 x d+1)
    l_M = inv(L_M)*V*v; %dx1            % INVERSE HERE WILL SLOW US DOWN!
    L_M(end+1,:) = l_M'; % new row 1xd, L_M now (d+1 x d)
    L_M(end,end+1) = sqrt(sigmasq + v'*v - l_M'*l_M); % bottom rt corner scalar
    l_Mi = L_M(end,end);

    temp = (inv(L_M)*V_prime);
    w = temp(end,:)';
    q = q + w.^2;


    % UPDATE Beta and MU
    Beta(end+1,1) = 1/l_Mi * (v'*Y - Beta'*l_M);

    mu = V_prime'*inv(L_M)'*Beta; % do this every loop? ANOTHER INVERSE!?
    % potentially a simpler form in the Appendix with m'?

    %%%%% Rename Variables
    V = V_prime;
    clear V_prime
    
    % UPDATE M (not needed every loop)
    %M = sigmasq*eye(length(I)) + V*V'; %not in the appendix, general eqn


end
% end the while loop, should have a final active set I


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOW DO FULL GPR

%%%%%% Execute FULL GPR

% FULL GPR
% Use either cholesky decomp -->    calcGPR_chol()
% or regular inverse with inv() --> calcGPR()
GPR_full = calcGPR(X,Y,nnum,"inv");
%ylim([-4, 4]);

% FULL GPR ON DOWNSAMPLED DATA
% Create a Downsampled Dataset N -> M
b = 4; % keep every jth point
[X_d,Y_d] = rawDownsample(X,Y,b);
GPR_d = calcGPR(X_d,Y_d,nnum,"inv");
%ylim([-4, 4]);



%%%%%%%%%%  CREATE SLIMDATA PLOTS

% Plot 1: Training DataPoints
 % Show pts in I (included)
 % Show pts in R (remaining)

f1 = figure;
plot(X(I),Y(I),'bo','markerfacecolor','b'), hold on
plot(X(R),Y(R),'ko'), hold off, grid on
xlabel('Training Data X Values'),ylabel('Training Data Y Values')
legend('Slim Data','Excess Raw Data')
title_str = strcat("Approx GPR on Slim Data (", ...
    num2str(length(X)),"raw to ", ...
    num2str(length(I)),"slim)");
title(title_str), grid on
%f.WindowState = 'maximized'; %make it full screen


% PREDICTION Points!
f2 = figure;
plot(X(I),Y(I),'bo','markerfacecolor','b'), hold on
xlabel('Training Data X Values'),ylabel('Training Data Y Values')

% k prediction points
X_beg = -nnum; X_end = nnum; 
X_star = [(-15+X_beg):2:(15+X_end)];

for k = 1:length(X_star)

    l_star = inv(L)*(K_Function(X_star(k),X(I)))';
    l_Mstar = inv(L_M)*l_star;

    Y_mu_star(k) = l_Mstar'*Beta;
    Y_sigmasq_star(k) = K_Function(X_star(k), X_star(k)) + ...
        - norm(l_star)^2 + ...
        sigmasq*norm(l_Mstar)^2;
end

% Plot of GPR w/ SlimData
errorbar(X_star,Y_mu_star,Y_sigmasq_star,'ro')
legend('Slim Data','Prediction \mu and \sigma^2')
title_str = strcat("Approx GPR on Slim Data (", ...
    num2str(length(X)),"raw to ", ...
    num2str(length(I)),"slim)");
title(title_str), grid on




%%%%%%%%%%%%%%%%%%%%%%%% METRICS OF FIT!
% Comparing SlimData to Downsampled Data

% Downsampled vs. Full 
Metrics.Downsample.SSD = sum((GPR_full.Y_Star_Hat - GPR_d.Y_Star_Hat).^2);
Metrics.Downsample.SST = sum((GPR_full.Y_Star_Hat - ones(GPR_full.npredict,1)*mean(GPR_full.Y_Star_Hat)).^2);
Metrics.Downsample.Rsqr = 1 - Metrics.Downsample.SSD / Metrics.Downsample.SST;

% SlimData vs. Full
Metrics.Slimdata.SSD = sum((GPR_full.Y_Star_Hat - Y_mu_star').^2);
Metrics.Slimdata.SST = sum((GPR_full.Y_Star_Hat - ones(GPR_full.npredict,1)*mean(GPR_full.Y_Star_Hat)).^2);
Metrics.Slimdata.Rsqr = 1 - Metrics.Slimdata.SSD / Metrics.Slimdata.SST;

disp(Metrics.Downsample.Rsqr);
disp(Metrics.Slimdata.Rsqr);
