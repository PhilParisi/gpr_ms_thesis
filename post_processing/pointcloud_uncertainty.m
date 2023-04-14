% Pointcloud Uncertainty Exploration
% taking a look at the uncertainty values across a given method

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


% basic stats about the pointcloud
exact_name = "exact_100_inference";
exact_num_pcds = length(dir(fullfile(strcat("..\results\PCD\",exact_name,"\block_size_400\predictions"), '*.pcd')));
exact_pcd = [];

% go thru the exact PCDs and combine into one
for i = 0:(exact_num_pcds-1)
    exact_cloud_path = strcat("..\results\PCD\",exact_name,"\block_size_400\predictions\cloud_",num2str(i),".pcd");
    exact_temp = double(loadpcd(exact_cloud_path))';
    exact_pcd = [exact_pcd; exact_temp];
end

%% plot the uncertainties


% Define the colormap you want to use
cmap = jet(256); % for example, using the jet colormap with 256 levels

% Plot the scatter points
scatter(exact_pcd(:,1), ...
        exact_pcd(:,2), ...
        0.1, ...
        exact_pcd(:,4))

% Add a colorbar to the plot
colormap(cmap);
c = colorbar;
c.Label.String = 'Uncertainty';

% Adjust the figure layout
axis equal
xlabel('x')
ylabel('y')

%%



% what are we comparing today? names of folders
exact_name = "exact_100_inference";
exact_num_pcds = length(dir(fullfile(strcat("../results/PCD/",exact_name,"/block_size_400/predictions"), '*.pcd')));

approx_name = "kmeans_90_inference_better"; % CHANGE THIS TO THE NAME OF THE FOLDER WTIH YOUR PCDS TO COMPARE w/ EXACT!
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


