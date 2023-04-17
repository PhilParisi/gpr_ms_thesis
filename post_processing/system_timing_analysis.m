% System Timing Analysis
    % go thru the txt files in results 
    % extract the times (time to run entire MPGPR)
    % plot graphs for comparison

    % the num_pcds (num_tiles) is from training_hp.xlsx
    % the .txts are from running the bash script gpr_automate.sh that runs the gpr bag mapper executables

clc, clearvars, close all

% add paths
if ispc()
    addpath("..\results\SysTiming\")
    dir_name = '..\results\SysTiming\';
else
    addpath("../results/SysTiming/")
    dir_name = '../results/SysTiming/';
end


file_list = dir([dir_name '*.txt']);


% Initialize the table to store the values and filenames
data_table = table('Size', [length(file_list), 2], 'VariableTypes', {'double', 'string'}, 'VariableNames', {'Value', 'Filename'});

% Loop over each file and extract the value
for i = 1:length(file_list)
    % Load the file
    file_data = importdata([dir_name file_list(i).name]);
    
    % Store the value and filename in the table
    data_table(i,:) = {file_data, file_list(i).name};
end

% store into individual arrays
average.times = data_table.Value(1:5)';
average.s = [10, 20, 25, 33, 50];
average.num_tiles = [332,323,314,305,230];
average.times_norm = average.times ./ average.num_tiles;

hybrid.times = data_table.Value(6:14)';
hybrid.s = 10:10:90;
hybrid.num_tiles = [344,346,340,342,335,323,314,312,222];
hybrid.times_norm = hybrid.times ./ hybrid.num_tiles;

kmeans.times = data_table.Value(15:23)';
kmeans.s = 10:10:90;
kmeans.num_tiles = [348,344,344,342,339,330,317,317,234];
kmeans.times_norm = kmeans.times ./ kmeans.num_tiles;

neighbor.times = data_table.Value(24:32)';
neighbor.s = 10:10:90;
neighbor.num_tiles = [367,343,350,318,343,302,307,166,167];
neighbor.times_norm = neighbor.times ./ neighbor.num_tiles;

random.times = data_table.Value(33:41)';
random.s = 10:10:90;
random.num_tiles = [353,344,344,340,339,331,321,303,236];
random.times_norm = random.times ./ random.num_tiles;

systematic.times = data_table.Value(42:46)';
systematic.s = [10, 20, 25, 33, 50];
systematic.num_tiles = [339,327,314,307,253];
systematic.times_norm = systematic.times ./ systematic.num_tiles;

exact.times = 186.043;
exact.s = 100;
exact.num_tiles = 350;
exact.times_norm = exact.times ./ exact.num_tiles;

disp('data loaded. run next section for plots')

%%
% Create Plot of Total Timing 
    % not normalized by number of processing tiles
clc, close all

figure()
plot(random.s, random.times, 'cs-'), hold on
plot(systematic.s, systematic.times, 'rv-')
plot(hybrid.s, hybrid.times, 'g^-')
plot(average.s, average.times, 'bx-')
plot(neighbor.s, neighbor.times, 'm*-')
plot(kmeans.s, kmeans.times, 'color',[1 0.5 0], 'marker', '+', 'linestyle','-')
plot(exact.s, exact.times,'ko','markerfacecolor','k')

legend('random','systematic','hybrid','average','neighbor','kmeans','exact')
set(gca,'xdir','reverse')
xlabel('downsample percent s')
ylabel('total system runtime (secs)')
grid on



% Create Plot of Normalized System Time (divide by num of pcd aka tiles generated)
figure()
plot(random.s, random.times_norm, 'cs-'), hold on
plot(systematic.s, systematic.times_norm, 'rv-')
plot(hybrid.s, hybrid.times_norm, 'g^-')
plot(average.s, average.times_norm, 'bx-')
plot(neighbor.s, neighbor.times_norm, 'm*-')
plot(kmeans.s, kmeans.times_norm, 'color',[1 0.5 0], 'marker', '+', 'linestyle','-')
plot(exact.s, exact.times_norm,'ko','markerfacecolor','k')


legend('random','systematic','hybrid','average','neighbor','kmeans','exact')
set(gca,'xdir','reverse')
xlabel('downsample percent s')
ylabel('normalized system runtime (secs)')
grid on