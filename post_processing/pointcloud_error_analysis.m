% Pointcloud Error Analysis
% for comparing an exact method with all data points to an approximate
% method with downsampled data points

% ros script outputs a bunch of pcd's, need to assemble those and compare
% only the pcd's that are held in common

% currently only works on Ubuntu sorry players

clc, clearvars, close all

%%%%%% SETUP

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

approx_name = "neighbor_10_inference";
approx_num_pcds = length(dir(fullfile(strcat("../results/PCD/",approx_name,"/block_size_400/predictions"), '*.pcd')));


%%%%%% SORT THE PCDs TO ISOLATE ONLY THE TILES IN COMMON 

matches = 0;
exact_pcd = [];
approx_pcd = [];

% go thru the approx PCD
for i = 0:(approx_num_pcds-1)


    % load approx PCD
    approx_cloud_path = strcat("../results/PCD/",approx_name,"/block_size_400/predictions/cloud_",num2str(i),".pcd");
    approx_temp = double(loadpcd(approx_cloud_path))';
    approx_temp_xy = approx_temp(:,1:2);


    % find the matching exact PCD
    for j = 0:(exact_num_pcds-1)

    % load exact PCD
    exact_cloud_path = strcat("../results/PCD/",exact_name,"/block_size_400/predictions/cloud_",num2str(j),".pcd");
    exact_temp = double(loadpcd(exact_cloud_path))';
    exact_temp_xy = exact_temp(:,1:2);


    % compare temp PCDs 
    testing = approx_temp_xy == exact_temp_xy;

    if sum(sum(testing)) == 2*length(testing)
        % we have a match!

        % add point cloud to their commanders
        exact_pcd = [exact_pcd; exact_temp];
        approx_pcd = [approx_pcd; approx_temp];
        matches = matches + 1;
        break % get out of searching for the exact pcd

    else
        % no match :{
        % keep looping through exact PCDs
        
    end

    end

end

disp(strcat("done matching clouds! match #: ", num2str(matches)))


%%%%%% CALCULATE DIFFERENCES

% Root Mean Square Error
diff = approx_pcd(:,3) - exact_pcd(:,3);
sqdiff = diff.^2;

RMSE = sqrt( mean(sqdiff) );


% Mean Absolute error
abs_diff = abs(approx_pcd(:,3) - exact_pcd(:,3));
MAE = mean(abs_diff,'all');


