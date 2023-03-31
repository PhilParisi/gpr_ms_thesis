% The Script that Makes all the Error Plots
% the data below is all from training_hp.xlsx

clc, clearvars, close all

%%%%%% EXACT
exact.dsample = 100;
exact.rawtime = 92.06;
exact.tiletime = 0.263;
exact.RMSE = 0;


%%%%%% AVERAGE DOWNSAMPLE
average.dsample = [50, 33, 25, 20, 10];
average.RMSE = [0.0580, 0.0683, 0.0847, 0.1047, 0.1663];
average.MAE = [0.0355, 0.0412, 0.0460, 0.0531, 0.0953];
average.rawtime = [44.53, 28.54, 22.69, 18.71, 10.74]; 
average.tiletime = [0.134,0.088,0.072,0.061,0.047];


%%%%%% HYBRID DOWNSAMPLE
hybrid.dsample = [90, 80, 70, 60, 50, 40, 30, 20, 10];
hybrid.RMSE = [0.0851, 0.1047, 0.0727, 0.0851, 0.0979, 0.0964, 0.0864, 0.1162, 0.1659];
hybrid.MAE = [0.0522, 0.0595, 0.0436, 0.0509, 0.0576, 0.0586, 0.0526, 0.0676, 0.1009];
hybrid.rawtime = [110.15, 89.53, 73.28, 65.94, 54.77, 43.44, 30.46, 20.50, 10.83]; 
hybrid.tiletime = [0.372, 0.309, 0.245, 0.227, 0.202, 0.161, 0.110, 0.083, 0.056];


%%%%%% KMEANS DOWNSAMPLE
kmeans.dsample = [90:-10:10];
kmeans.RMSE = [0.1281, 0.1085,0.1127,0.1243,0.1111,0.0949,0.1060,0.1010,0.1513];
kmeans.MAE = [0.0729,0.0619,0.0632,0.0709,0.0635,0.0565,0.0624,0.0626,0.0939];
kmeans.rawtime = [150.76,125.29,103.57,82.77,67.30,49.68,35.69,21.34,10.29]; 
kmeans.tiletime = [0.857,0.712,0.606,0.473,0.391,0.305,0.215,0.133,0.074];


%%%%%% NEIGHBOR DOWNSAMPLE
neighbor.dsample = [90:-10:10];
neighbor.RMSE = [0.1213,0.1308,0.1219,0.1636,0.1250,0.1348,0.1402,0.1631,0.1979];
neighbor.MAE = [0.0723,0.0814,0.0749,0.0919,0.0728,0.0828,0.0917,0.1094,0.1376];
neighbor.rawtime = [158.00,119.32,97.72,72.81,56.88,39.44,27.93,17.03,9.90]; 
neighbor.tiletime = [0.836,0.645,0.534,0.394,0.320,0.219,0.154,0.101,0.061];


%%%%%% SYSTEMATIC DOWNSAMPLE
systematic.dsample = [50, 33, 25, 20, 10];
systematic.RMSE = [0.0564,0.0674,0.0718,0.0874,0.2134];
systematic.MAE = [0.0362,0.0393,0.0422,0.0512,0.1034];
systematic.rawtime = [43.43,30.45,23.47,19.18,11.16]; 
systematic.tiletime = [0.128,0.093,0.075,0.063,0.044];


%%%%%% RANDOM DOWNSAMPLE
random.dsample = [90:-10:10];
random.RMSE = [0.0517,0.0543,0.0744,0.0665,0.0638,0.0625,0.0830,0.0830,0.1819];
random.MAE = [0.0319,0.0346,0.0408,0.0370,0.0385,0.0384,0.0446,0.0477,0.0849];
random.rawtime = [84.08,72.50,63.67,53.40,45.43,36.81,28.65,18.99,11.00]; 
random.tiletime = [0.238,0.211,0.185,0.157,0.134,0.111,0.089,0.063,0.047];



% mega data
mega.average = average;
mega.exact = exact;
mega.hybrid = hybrid;
mega.kmeans = kmeans;
mega.neighbor = neighbor;
mega.random = random;
mega.systematic = systematic;

disp('results loaded')


%% PLOTS
clc, close all

% Accuracy vs. Time
%figure()



% Accuracy vs. Downsample
figure()
%for i = 1:6
    %plot(exact.dsample,exact.RMSE,'ko','markerfacecolor','k'), hold on
    plot(average.dsample,average.RMSE,'bx-'), hold on
    plot(systematic.dsample,systematic.RMSE,'rv-')
    plot(hybrid.dsample,hybrid.RMSE,'g^-')
    plot(neighbor.dsample,neighbor.RMSE,'m*-')
    plot(random.dsample,random.RMSE,'cs-')
    plot(kmeans.dsample,kmeans.RMSE, 'color', [1 0.5 0], 'marker', '+')
%end

title('Error (RMSE) vs. Downsample')
xlabel('Downsample Percent \it{s}')
ylabel('Root Mean Square Error (m)')
grid on
xlim([0 100])
legend('average','systematic','hybrid','neighbor','random','kmeans','location','northwest')
set(gca,'XDir','reverse')
hold off

% Time vs. Downsample

figure()
%for i = 1:6
    plot(exact.dsample,exact.tiletime,'ko','markerfacecolor','k'), hold on 
    plot(average.dsample,average.tiletime,'bx-')
    plot(systematic.dsample,systematic.tiletime,'rv-')
    plot(hybrid.dsample,hybrid.tiletime,'g^-')
    plot(neighbor.dsample,neighbor.tiletime,'m*-')
    plot(random.dsample,random.tiletime,'cs-')
    plot(kmeans.dsample,kmeans.tiletime, 'color', [1 0.5 0], 'marker', '+')
%end

title('Mean Processing Time vs. Downsample')
xlabel('Downsample Percent \it{s}')
ylabel('Mean Processing Time (secs)')
xlim([0 100])
grid on
legend('exact','average','systematic','hybrid','neighbor','random','kmeans')
set(gca,'XDir','reverse')
hold off