%%%%%% Datasets for GPR

% Purpose: multiple sections that can generate new random datasets
    % for training on w/ GPR

%% 2D LARGE SEAMOUNTS
clc, clear all, close all

% this dataset is designed to include a flat seafloor
% with two salient features
% Goal is to test how the downsampling algorithms sample the features!

num_pts = 20000; %takes about 20 seconds to do inverse

data.x = 0:20000;
data.length = length(data.x);

% First Quarter, flat seafloor
pts1 = 3000;
data.y = ones(1,data.length);
data.y(1:pts1-1) = data.y(1:pts1-1) + 0.1*sin(data.x(1:pts1-1)/data.length*100);

% Flattop Seamount
pts2 = 7000;
data.y(pts1:pts2) = data.y(pts1:pts2) + 5;

% Gaussian Seamount
pts3 = 7000; pts4 = 19000;
mu = (pts4+pts3)/2; sigma = (pts4-pts3)/8;
data.y(pts3:pts4) = 10000*1/(sqrt(2*pi) * sigma) * exp(-0.5* (1/sigma*(data.x(pts3:pts4)-mu)).^2 ) + ...
    data.y(pts3:pts4);

% Overall Noise (sensor error, gaussian)
noise = randn(1,num_pts+1)/70;
data.y = data.y + noise;

% Make Figure
figure
plot(data.x,data.y,'.','MarkerSize',1)
ylim([0 10])
xlabel('Index'),ylabel('Depth')
title_str = strcat("Raw Data with ",num2str(data.length), " points");
title(title_str)

%%%%%%%%%%%%%

% Rename Key variables and chop off last value (gives a nice even number)
X = data.x(1:end-1); Y = data.y(1:end-1);

% Make them VERTICAL (nx2)
X = X'; Y = Y';

% Save workspace variables to .mat
filename = strcat("GPR_dataset_Seamounts_",...
    num2str(length(X)),...
    "pts_2D.mat");
save(filename,"X","Y");

%% 2.5D SAND RIPPLES
clc, clear all, close all

% dimensions of map
map_size = 50; % keep this even if you can 
x = (-map_size/2:(map_size/2-1)); % + randn(1,map_size)/2;
y = (-map_size/2:(map_size/2-1)); % + randn(1,map_size)/2;
[X,Y] = meshgrid(x,y);

% 2Matrices map_size x map_size
X = X + randn(map_size,map_size);
Y = Y + randn(map_size,map_size);
Z = 0.5*sin(X) + Y/20;
Z = Z + randn(map_size,map_size)/5;

% Reshape to Plot (map_size^2 x 1)
X1 = reshape(X,map_size^2,1);
Y1 = reshape(Y,map_size^2,1);
Z1 = reshape(Z,map_size^2,1);

% Plot
scatter3(X1,Y1,Z1,'.')
xlabel('X'),ylabel('Y'),zlabel('Depth')
zlim([-5 5])


% Save workspace variables to .mat
filename = strcat("GPR_dataset_SandRipples_",...
    num2str(length(X1)),...
    "pts_2p5D.mat");
save(filename,"X1","Y1","Z1");


%% INVERSE CALCULATION SPEED DIAGRAM
clc, close all, clear all

% How long does it take matlab to do an inverse?

pts = [1000 5000 10000 15000 20000];
time_calc = zeros(1,length(pts));
for i = pts
    jj = rand(pts(i),pts(i));
    tic
    anss = inv(jj);
    time_calc(i) = toc
end

%%
close all
figure
plot(pts,time_calc,'--b'), hold on
plot(pts,time_calc,'ob','MarkerFaceColor','b')
title('Time to Calculate Inverse of a Matrix')
ylabel('Calc Time (secs)'),xlabel('# Pts in Matrix')
grid on
xticks = pts; xticks(1) = 0;
set(gca,'XTickLabel',xticks)
