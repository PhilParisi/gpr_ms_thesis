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

%data.Z = ones(height(data),1)*100;

% parameters
range_var = .10 / 6; % 10cm rating at 6 sigma = 0.0167m
        % or is it +/- 10cm rating = 20cm. 20cm at 6 sigma = 0.033m
nbeams = 256; swath_angle = 120*pi/180; % CHECK THESE VALUES
angle_var = swath_angle / nbeams;

% calculate variances
ranges = sqrt(data.Y.^2 + data.Z.^2);
angles = atan(abs(data.Y)./abs(data.Z));

variances = (cos(angles).*range_var).^2 + (ranges.*sin(angles).*angle_var).^2;


% constant uncertainties
const_variances = ones(length(variances))*0.04;

cupper = data.Z + 2*sqrt(const_variances);
clower = data.Z - 2*sqrt(const_variances);


% PLOTS
figure()
%subplot(1,2,1)
plot(data.Y,data.Z, 'k.'), hold on

upper = data.Z + 2*sqrt(variances);
lower = data.Z - 2*sqrt(variances);

%plot uncerts
plot(data.Y,upper,'b-.')

plot(data.Y,cupper,'r-')
plot(data.Y,clower,'r-')

plot(data.Y,lower,'b-.')

xlim([-10 10])
grid on
title('95% Confidence [+/-2\sigma] Intervals on a Sonar Ping'),
xlabel('Position (m)'), ylabel('Depth (m)')



legend('raw data','Changing Uncertainty','Constant Uncertainty','location','northwest')

%{
% constant uncertainties
const_variances = ones(length(variances))*0.04;
%subplot(1,2,2)
plot(data.Y,data.Z, 'k.'), hold on

upper = data.Z + 2*sqrt(const_variances);
lower = data.Z - 2*sqrt(const_variances);
plot(data.Y,upper,'b-')
plot(data.Y,lower,'b-')
xlim([-10 10])
grid on
title('2*std-dev bounds with constant uncertainties over a ping'),
xlabel('Position (m)'), ylabel('Depth (m)')
legend('raw data','2std bounds')
%}