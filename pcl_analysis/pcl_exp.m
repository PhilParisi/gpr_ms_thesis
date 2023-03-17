% PCL Analysis 

if isunix()
    addpath("../data/")
    addpath("../pcl_functions/")
else
    addpath("..\data\")
    addpath("..\pcl_functions\")
end

clc, clearvars, close all, format compact


% load a pointcloud
cloud = loadpcd('cloud_1.pcd');
cloud = double(cloud);

% the likelihood function
pt1 = cloud(:,1); % approx
pt2 = cloud(:,2); % exact
    
    % note that we already have sigma^2 (variances)
likelihood_num = exp(-0.5 * (pt1(3) - pt2(3))^2 / (sqrt(pt1(4)) + sqrt(pt2(4))) );
likelihood_den = sqrt(2*pi * (sqrt(pt1(4)) + sqrt(pt2(4)))  );
ptLikelihood = likelihood_num / likelihood_den


likelihood_num = exp(-0.5 * (7.5 - 8)^2 / (0.1 + 0.1) );
likelihood_den = sqrt(2*pi * (0.1 + 0.1) );
ptLikelihood = likelihood_num / likelihood_den
