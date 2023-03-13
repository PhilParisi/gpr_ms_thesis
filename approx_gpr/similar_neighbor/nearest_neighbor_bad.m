% Nearest Neighbor
% main assumption: we want points that are further apart
% calculate distances, select largest ones

clc, clearvars, close all

% Raw Data
if ispc()
    addpath("..\..\data\")
else
    addpath("../../data/")
end

rawdata = readtable("wiggles_single_ping.csv");
data = rawdata(77:end-37,[2 3]);

% Params
dsample_percent = 0.3;
sample_vals = dsample_percent * height(data);


% Loop over every datapoint


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

    D(i) = mean(distances)/2;

end

[D2, idx] = sort(D,'descend');



plot(data.Y,data.Z,'k.'), hold on
plot( data.Y(idx(1:sample_vals)), data.Z(idx(1:sample_vals)),'ro' )

%%

clc, clearvars, close all

probs = [0.5, 0.3, 0.1, 0.05, 0.05];

% Set the weights for the randsample function
%weights = [0.5, (1 - 0.5) * probs(2:end) / sum(probs(2:end))];
weights = probs;

% Randomly select from the probs array using the weights
selected_value = randsample(length(probs), 2, true, weights);

% Display the selected value
disp(probs(selected_value));