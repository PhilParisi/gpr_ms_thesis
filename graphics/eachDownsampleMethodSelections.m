% The purpose of this script to to a create a graphic
% that shows the selection points for each downsample method

% Random, Systematic, Hybrid, Point Average, Info Gain
% We're going to do a 10% downsample so it's super clear

% Please know this code is not super optimized so. Just for a few figures.

clc, clearvars, close all

% Raw Data
if ispc()
    addpath("..\data")
else
    addpath("../data")
end

rawdata = readtable("wiggles_single_ping.csv");
data = rawdata(77:end-37,[2 3]);

% Params
dsample_percent = 0.2; % 0.1, 0.2, 0.5 are best, others have to do rounding so each method may have a different number of points
sample_vals = round(dsample_percent * height(data));


%%%%%%%%%% DOWNSAMPLE METHODS

% systematic
j = round(1/dsample_percent);
sys_start = 3;
sys_indices = sys_start:j:height(data);

% random
perm = randperm(height(data)); % random permutation
rand_indices = perm(1:sample_vals);

% hybrid
h = 2*(1/dsample_percent);
hyb_start = 17;
hyb_1 = hyb_start:h:height(data); % systematic
hyb_perm = randperm(height(data));
hyb_2 = hyb_perm(1:(sample_vals-length(hyb_perm)));
hyb_ind = [hyb_1 hyb_2];

while length(unique(hyb_ind)) ~= sample_vals
    hyb_perm = randperm(height(data));
    hyb_2 = hyb_perm(1:(sample_vals-length(hyb_1)));
    hyb_ind = [hyb_1 hyb_2];
end
%hyb_2 = [92 23 59 83 24]; %got these thru random generation in above 2 lines

hyb_indices = hyb_ind;


% averaging
p = round(1/dsample_percent);
counter = 1;
for i = 1:p:height(data)        % loop over all points
    if (i + p - 1) <= height(data)
        
        avg_data.Y(counter) = sum(data.Y(i:i+p-1)) / p;
        avg_data.Z(counter) = sum(data.Z(i:i+p-1)) / p;
        counter = counter + 1;
    end
end


% dissimilar neighbor
% go over every point
for i = 1:height(data)

    % indices to check
    if i == 1
        explore = [i+1, i+2];
    elseif i == height(data)
        explore = [i-1, i-2];
    else
        explore = [i-1, i+1];
    end

    % calc distances to explore indices
    distances(1) = sqrt( (data.Y(explore(1)) - data.Y(i))^2   + (data.Z(explore(1)) - data.Z(i))^2 );
    distances(2) = sqrt( (data.Y(explore(2)) - data.Y(i))^2   + (data.Z(explore(2)) - data.Z(i))^2 );
    distances = sort(distances); %smallest first

    % calc infogain
    %infogain(i) = log(distances(2) / distances(1));
    infogain(i) = (distances(1) / distances(2)); %small over big

end


% sort the info gain
[sorted_infogain, idx] = sort(infogain, 'ascend'); %smaller first (we want those)

infogain_pts = data(idx(1:round(sample_vals)),["Y", "Z"]);


% kmeans
kdata = table2array(data);
[idx,C] = kmeans(kdata,sample_vals);



%{
% ONE PLOT
zshift = 0.1;

figure()

% raw data
plot(data.Y, data.Z, 'k.'), hold on

% systematic
plot(data.Y(sys_indices), data.Z(sys_indices) + zshift, 'rs')

% random
plot(data.Y(rand_indices), data.Z(rand_indices) + 2*zshift, 'go')

% hybrid
plot(data.Y(hyb_indices), data.Z(hyb_indices) + 3*zshift, 'md')

% average
plot(avg_data.Y, avg_data.Z + 4*zshift, 'ch')

% infogain
plot(infogain_pts.Y, infogain_pts.Z + 5*zshift, 'bv')

% kmeans
plot(C(:,1),C(:,2),'gx')

% extras
title('Selection Points')
xlabel('Position (m)')
ylabel('Depth (m)')
grid on
legend('raw data','systematic', 'random', 'hybrid', 'average', 'info gain',...
    'location','northwest')
%}


% SUBPLOTS

figure()

% raw + systematic data
subplot(3,2,2)
plot(data.Y, data.Z, 'k.'), hold on
plot(data.Y(sys_indices), data.Z(sys_indices), 'rs','MarkerFaceColor','r'), hold off
title('Systematic Inclusions'), xlabel('Position (m)'), ylabel('Depth (m)')
grid on

% raw + random
subplot(3,2,1)
plot(data.Y, data.Z, 'k.'), hold on
plot(data.Y(rand_indices), data.Z(rand_indices), 'o','Color',[0 0.5 0],'MarkerFaceColor',[0 0.5 0]), hold off
title('Uniform Random Inclusions'), xlabel('Position (m)'), ylabel('Depth (m)')
grid on

% raw + hybrid
subplot(3,2,3)
plot(data.Y, data.Z, 'k.'), hold on
plot(data.Y(hyb_indices), data.Z(hyb_indices), 'md','MarkerFaceColor','m'), hold off
title('Hybrid Inclusions'), xlabel('Position (m)'), ylabel('Depth (m)')
grid on

% raw + average
subplot(3,2,4)
plot(data.Y, data.Z, 'k.'), hold on
plot(avg_data.Y, avg_data.Z, '^', 'Color', [0.5 0 0.5], 'MarkerFaceColor', [0.5 0 0.5]), hold off
title('Point Averaging Inclusions'), xlabel('Position (m)'), ylabel('Depth (m)')
grid on

% raw and infogain
subplot(3,2,5)
plot(data.Y, data.Z, 'k.'), hold on
plot(infogain_pts.Y, infogain_pts.Z,'bv','MarkerFaceColor','b')
title('Dissimilar Neighbor Inclusions'), xlabel('Position (m)'), ylabel('Depth (m)')
grid on

% raw and kmeans
subplot(3,2,6)
plot(data.Y, data.Z, 'k.'), hold on
plot(C(:,1),C(:,2),'>','Color',[1 0.5 0], 'MarkerFaceColor', [1 0.5 0])
title('K-means Inclusions'), xlabel('Position (m)'), ylabel('Depth (m)')
grid on
