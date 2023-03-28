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


%%

% add paths
if ispc()
    addpath("..\results\PCD\")
    addpath("..\pcl_functions\")
else
    addpath("../results/PCD/")
    addpath("../pcl_functions/")
end

% what are we comparing today? names of folders
exact_name = "exact_100_inference";
exact_num_pcds = length(dir(fullfile(strcat("../results/PCD/",exact_name,"/block_size_400/predictions"), '*.pcd')));

approx_name = "average_50_inference";
approx_num_pcds = length(dir(fullfile(strcat("../results/PCD/",approx_name,"/block_size_400/predictions"), '*.pcd')));

%%
exact_pcd = [];

% load exact PCDs
for i = 0:(exact_num_pcds-1)
    exact_cloud_path = strcat("../results/PCD/",exact_name,"/block_size_400/predictions/cloud_",num2str(i),".pcd");
    temp_pcd = double(loadpcd(exact_cloud_path));
    exact_pcd = [exact_pcd; temp_pcd'];
end


approx_pcd = [];
% load approx PCDs
for i = 0:(approx_num_pcds-1)
    approx_cloud_path = strcat("../results/PCD/",approx_name,"/block_size_400/predictions/cloud_",num2str(i),".pcd");
    temp_pcd = double(loadpcd(approx_cloud_path));
    approx_pcd = [approx_pcd; temp_pcd'];
end