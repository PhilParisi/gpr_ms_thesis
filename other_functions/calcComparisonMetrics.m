function [Metrics] = calcComparisonMetrics(GPR,GPR_d)
% calcComparisonMetrics compares two GPR results (fit of prediction means)
% Input Full GPR solution first, Slim GPR second

Mean1 = GPR.Y_Star_Hat; Mean2 = GPR_d.Y_Star_Hat;

if (length(Mean1) ~= length(Mean2))
    disp('Means are not of equal length. They must be!')
end

%%%%%% Calculate Metrics
% each method against TRAINING DATA? GRIDDED PRODUCT?

%%% Sum of Squared Differences - between exact and approx
Metrics.SSD = sum((Mean1 - Mean2).^2);

%%% R-Squared - between exact and approx
Metrics.SST = sum((Mean1 - ones(length(Mean1),1)*mean(Mean1)).^2);
Metrics.Rsqr = 1 - Metrics.SSD / Metrics.SST;

%%% Calculate Speed Comparison
%Metrics.PerFaster = GPR.compute_time / GPR_d.compute_time;

%%% Output Results to Command Window
fprintf("RSqr (Exact-Approx): %0.3f\n",Metrics.Rsqr)
%fprintf("Times Faster (Exact-Approx): %0.3f\n",Metrics.PerFaster)

%%% Calculate Speed Comparison
%Metrics.PerFaster = GPR.compute_time / GPR_d.compute_time;


% Plots!
calcComparisonPlots(GPR,GPR_d);


end

