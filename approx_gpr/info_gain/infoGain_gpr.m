% Fast Forward Selection to Speed up Sparse GPR
% Authors: Seeger, Williams, Lawrence
% Programmer: Parisi
% Date: 18May2022

% GOAL
    % Input: raw training data of n total points (inputs and outputs)
    % Output: down-selected subset of d points from n total pts
        % these will then be used as the inputs for the GPR model

% This script should perform the repeated selection of slim data points
% information gain will be calculated, new point added to slim data, repeat
% model updates are needed, though this is not fully understood yet
% are updates done to hyperparameters? the likelihood assumptions? idk tbh


%%
clc
if ispc()
    addpath("..\..\gpr_functions\")
else
    addpath("../../gpr_functions/")
end



%%
%%%%%%%%% Bring in Raw Data
clc, clear all, close all, format compact, %close all

%load("../../data/approx_methods_datasets/matlab_gpr_datasets/A_random_training_data_gpr_40.mat")
%load("~/Data/approx_methods_datasets/matlab_gpr_datasets/B_random_training_data_gpr_5000.mat")
%nnum = length(X);
%addpath("..\..\data")
%load("..\..\data\wiggles_sample_ping.csv")
data = readtable("wiggles_single_ping.csv");

%disp(strcat("...training data with ", num2str(nnum), " datapoints loaded..."))
%clearvars -except nnum X Y

X = data.X; Y = data.Z; nnum = length(X);

figure(1), plot(X,Y,'b.','markersize',5)

%% SELECT MULTIPLE POINTS (GROUP OF INCLUSION POINTS)

clc, clearvars -except nnum X Y, close all
%Downsampling as Described in Paper

% total points wanted in slim data [threshold for ending algo]
totalpts_slimdata = round(0.2*nnum);

%%%%%%%%%  Hyperparamters
hp.L = 10;
hp.sigma_p = 1.2;
hp.sigma_n = 0.5;
hp.kerneltype = 'sparse';       %exact or sparse

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
K_I_I = K_Function(X(I),X(I),hp);      % Self Covariance, I vs. I, dxd
K_I_dot = K_Function(X(I),X(:),hp);    % Joint Covariance, I vs. All, dxn
P_I = inv(K_I_I) * K_I_dot;         % !?!?INVERSE!?!? dxd * dxn

% cholesky decomp of K_I_I
L = chol(K_I_I,'lower');        % dxd    % Lower Triangular
V = inv(L)*K_I_dot;             % dxn
%Ident_I = eye(length(I));       % dxd    % Identity Matrix
M = hp.sigma_p^2*eye(length(I)) + V*V';     % dxd    % Noisy Covariance maybe? 

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
        l_i = sqrt(round(K_Function(X(R_ind),X(R_ind),hp) - p_i, 5)); 

        epsilon_i = 1/((hp.sigma_p/l_i)^2 + 1 - q_i);
        kappa_i = epsilon_i*(1+2*(hp.sigma_p/l_i)^2);

        % info gain equation
        info_gain(R_ind,1) = -log(hp.sigma_p/l_i) + ... % NATURAL LOG?
            -0.5*(log(epsilon_i) + ...
            epsilon_i*(1-kappa_i)*(hp.sigma_p^-2)*((Y(R_ind)-mu(R_ind))^2) + ...
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
    L(end,end+1) = sqrt( K_Function(X(inclusion_ind),X(inclusion_ind),hp) - p(inclusion_ind) ); % bottom rt scalar
    
    li = L(end,end); % bottom rt corner of Lprime
    vi = V(:,inclusion_ind); %dx1 col vector
    

    % UPDATE V (dxn) to V_prime (d+1 x n, add new row)
    v = 1/li * (K_Function(X,X(inclusion_ind),hp) - V'*vi); %nx1 col vector
    V_prime = [V; v'];


    p = p + v.^2; % is this the right approach?

    % UPDATE L_M (dxd to d+1 x d+1)
    l_M = inv(L_M)*V*v; %dx1            % INVERSE HERE WILL SLOW US DOWN!
    L_M(end+1,:) = l_M'; % new row 1xd, L_M now (d+1 x d)
    L_M(end,end+1) = sqrt(hp.sigma_p^2 + v'*v - l_M'*l_M); % bottom rt corner scalar
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


%%%%%%%%%%  CREATE PLOTS
close all

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





%%

%%%%% PREDICTIONS!
f2 = figure;
plot(X(I),Y(I),'bo','markerfacecolor','b'), hold on
xlabel('Training Data X Values'),ylabel('Training Data Y Values')

% k prediction points
X_beg = -nnum; X_end = nnum; 
X_Star = [(-15+X_beg):.2:(15+X_end)];

for k = 1:length(X_Star)

    l_star = inv(L)*(K_Function(X_Star(k),X(I),hp))';
    l_Mstar = inv(L_M)*l_star;

    Y_mu_star(k) = l_Mstar'*Beta;
    Y_sigmasq_star(k) = K_Function(X_Star(k), X_Star(k), hp) + ...
        - norm(l_star)^2 + ...
        hp.sigma_p^2*norm(l_Mstar)^2;
end
% 
% errorbar(X_star,Y_mu_star,Y_sigmasq_star,'ro')
% legend('Slim Data','Prediction \mu and \sigma')
% title_str = strcat("Approx GPR on Slim Data (", ...
%     num2str(length(X)),"raw to ", ...
%     num2str(length(I)),"slim)");
% title(title_str), grid on

Y_Star_Hat = Y_mu_star;
Y_Star_Var = Y_sigmasq_star;

sortobj = [X_Star', Y_Star_Hat', Y_Star_Var'];
sortobj = sortrows(sortobj);
sorted.X_Star = sortobj(:,1); sorted.Y_Star_Hat = sortobj(:,2); sorted.Y_Star_Var = sortobj(:,3);


f3 = figure;
p1 = plot(sorted.X_Star,sorted.Y_Star_Hat + 2*sqrt(sorted.Y_Star_Var),'r','LineWidth',2); hold on %upper bound
plot(sorted.X_Star,sorted.Y_Star_Hat - 2*sqrt(sorted.Y_Star_Var),'r','LineWidth',2); % lower bound
p2 = plot(sorted.X_Star,sorted.Y_Star_Hat,'r--','Linewidth',2); % prediction means
p3 = plot(X(I),Y(I),'bo','markerfacecolor','b'); % pts in I
p4 = plot(X(R),Y(R),'ko'); % pts in R
xlabel('X Values'), ylabel('Y Values'), title('Gaussian Process Regression')
grid on, legend([p3 p4 p2 p1],"Slim Data","Excess Raw Pts","Prediction \mu","Prediction 2\sigma")

