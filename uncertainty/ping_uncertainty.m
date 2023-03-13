% Calculating uncertainties along a ping
% Using simplified uncertainty model

clc, clearvars, close all

% add paths
if ispc()
    addpath("..\data")
else
    addpath("../data")
end

% load and organize data
rawdata = readtable("wiggles_single_ping.csv");
data = rawdata(:,[1 2 3]);
data.X = zeros(size(data.X));
data.Y = data.Y - data.Y(107); % center around middle point

data.Z = ones(height(data),1);

% parameters
range_var = .20 / 6; % +/- 10cm rating = 20cm. 20cm at 6 sigma = 0.033m
nbeams = 256; swath_angle = 90*pi/180; % CHECK THESE VALUES
angle_var = swath_angle / nbeams;

% calculate variances
ranges = sqrt(data.Y.^2 + data.Z.^2);
angles = atan(abs(data.Y)./abs(data.Z));

variances = (cos(angles).*range_var).^2 + (ranges.*sin(angles).*angle_var).^2;



% PLOTS
figure()
plot(data.Y,data.Z, 'k.'), hold on

upper = data.Z + 2*sqrt(variances);
lower = data.Z - 2*sqrt(variances);
plot(data.Y,upper,'b-')
plot(data.Y,lower,'b-')