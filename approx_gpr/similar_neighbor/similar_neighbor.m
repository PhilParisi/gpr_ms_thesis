% Information Gain from Information Theory


clc, clearvars, close all

% Raw Data
if ispc()
    addpath("..\..\..\data\")
else
    addpath("../../../data/")
end

rawdata = readtable("wiggles_single_ping.csv");
data = rawdata(77:end-37,[2 3]);

% Params
dsample_percent = 0.25;

sample_vals = dsample_percent * height(data);


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


plot(data.Y, data.Z, 'k.'), hold on
plot(infogain_pts.Y, infogain_pts.Z,'ro')



