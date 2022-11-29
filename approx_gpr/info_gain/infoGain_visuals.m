% Fast Forward Selection to Speed up Sparse GPR
% Authors: Seeger, Williams, Lawrence
% Programmer: Parisi
% Date: 18May2022

% GOAL
    % Input: raw training data of n total points (inputs and outputs)
    % Output: down-selected subset of d points from n total pts
        % these will then be used as the inputs for the GPR model

% This script is a basic implementation of the aforementioned paper
% With given input raw data, a few indices are randomly generated
% These become our 'slim dataset'.
% We then search through the remaining indices and calculate info_gain
% The pt with the highest info_gain is added to the 'slim dataset'

% This script only does the INFO GAIN CALCULATE. Most basic building block.
% see other scripts for more advacned implementation.

%% Add GPR Functions to Path
addpath("gpr_functions/")

%%
%%%%%%%%% Bring in Raw Data
clc, close all, clear all, format compact

load("~/Data/approx_methods_datasets/matlab_gpr_datasets/A_random_training_data_gpr_40.mat")
%load("~/Data/approx_methods_datasets/matlab_gpr_datasets/B_random_training_data_gpr_5000.mat")
nnum = length(X);

disp(strcat("...training data with ", num2str(nnum), " datapoints loaded..."))
clearvars -except nnum X Y

%% SELECT ONE POINT (NEXT INCLUSION POINT)
clc, clearvars I R, close all
%Downsampling as Described in Paper

%%%%%%%%% Hyperparamters
sigma = 2;
sigmasq = sigma^2;

%%%%%%%%% Indexing Variables
I = []; % Active Set of Indices to Include (starts empty), w/ d values
    % add to this list with --> I(end+1) = value
R = 1:nnum; % Remaining Indices not in Active Set, I
    % remove from this list with --> R(index) = []

slimdatapts = 5;
% let's just get some values into I for now! 3 of them
add_index = randi(nnum,1,slimdatapts);        % 3 random indices from training data
while length(unique(add_index)) < slimdatapts
    add_index = randi(nnum,1,slimdatapts);    % makes sure the 3 indices are different
end

I = [I add_index];
R(add_index) = [];



% now we can begin to calculate terms needed for Information Gain 
% calculate the BEFORE we consider i terms
%%%%%%%%% CURRENT TERMS

% beginning terms
K_I_I = K_Function(X(I),X(I));      % Self Covariance, I vs. I
K_I_dot = K_Function(X(I),X(:));    % Joint Covariance, I vs. All
P_I = K_I_I * K_I_dot;

% cholesky decomp of K_I_I
L = chol(K_I_I,'lower');        % dxd    % Lower Triangular
V = inv(L)*K_I_dot;             % dxn
Ident_I = eye(length(I));       % dxd    % Identity Matrix
M = sigmasq*Ident_I + V*V';     % dxd    % Noisy Covariance maybe? 

% cholesky decomp of M
L_M = chol(M,'lower');          % dxd
Beta = inv(L_M)*V*Y;            % dx1

% p and q
p = diag(V'*V);                 % nx1
q = diag(V'*inv(M)*V);          % nx1

% mu
mu = V'*inv(L_M')*Beta;         % IS THIS RIGHT??? ^-T in paper

%%%%%%%%%%  INFORMATION GAIN
% go thru each i in R, and calculate info gain (just 1 i for now...)
info_gain = [];
for ind = R       % go thru all R indices
    
    % APPROXIMATE INFORMATION GAIN
    p_i = p(ind);
    q_i = q(ind);
    l_i = sqrt(K_Function(X(ind),X(ind)) - p_i);
    epsilon_i = 1/((sigma/l_i)^2 + 1 - q_i); 
    kappa_i = epsilon_i*(1+2*(sigma/l_i)^2);
    info_gain(ind,1) = -log(sigma/l_i) + ...
        -0.5*(log(epsilon_i) + ...
        epsilon_i*(1-kappa_i)*(sigma^-2)*((Y(ind)-mu(ind))^2) + ...
        -kappa_i + 2);
               % recall log is natural log, log10 is log

   
end

close all
% Create Plot of Pts
f = figure;
%figure
subplot(2,1,2)
plot(info_gain,'bo'), hold on
plot(add_index,zeros(1,slimdatapts),'ko','markerfacecolor','k'), hold on
[include_pt_val, include_pt_ind] = max(info_gain);
plot(include_pt_ind,include_pt_val,'bo','markerfacecolor','b')
xlabel('Index'), ylabel('Information Gain'), grid on
legend('DataPoints','Pts Included Already','New Pt to Include')
title('Information Gained by Including a Given Index into Slim Data')
hold off

subplot(2,1,1)
%figure
plot(X,Y,'bo'), hold on
plot(X(add_index),Y(add_index),'ko','markerfacecolor','k'), hold on
plot(X(include_pt_ind),Y(include_pt_ind),'bo','markerfacecolor','b'), hold on
xlabel('Training Data X Values'),ylabel('Training Data Y Values')
title('Training DataPoints')
grid on
legend('DataPoints','Pts Included Already','New Pt to Include')
f.WindowState = 'maximized'; %make it full screen

