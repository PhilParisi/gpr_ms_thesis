% Uncertainty calculations for a given ping
clc, clearvars, close all

% Parameters
r = 10;
dr = 0.2;
fdr = dr/r*100; % fractional uncertainty

theta = 45*pi/180;
dtheta = 1/2*pi/180; 
fdtheta = dtheta/theta*100; % fractional

%%% Uncertainty - Method 1
% d = rj, j = cos(0)
j = cos(theta);
dj = abs(sin(theta))*dtheta;
fdj = dj/j*100;

% combine r and j (product, so add squares of f)
d = r*cos(theta);
fdd = sqrt( (fdj/100)^2 + (fdr/100)^2 )
dd = d*fdd