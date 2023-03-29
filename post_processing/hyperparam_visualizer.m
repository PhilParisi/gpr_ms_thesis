% Hyperpameters Visualizer
% the csv_gpr_gui outputs a csv with all the hp's tested
% this script opens the csv and find the max LML

clc, clearvars, close all

if ispc()
    addpath('..\results\LML')
else
    addpath('../results/LML/')
end

%method = ["random","hybrid","kmeans","neighbor"];
%percents = 10:10:90;
method = ["kmeans"];
percents = ["40_big","30_big","20_big","10_big",10,20,30,40,50,60,70,80,90];


for j = 1:length(method)
    for i = percents

        dpercent = num2str(i);

        % load the csv
        filename = strcat(method(j),"_",dpercent,"_lml");
        data = readtable(strcat(filename,".csv"));
        %data = readtable("exact_100_lml.csv");
        % remove rows that have NAN LML
        %to_remove = isnan(data.lml)';
        %data = data(~to_remove,:);

        % make a scatter plot
        data = sortrows(data,'lml','descend'); 
        final_ind = round(0.98*length(data.lml));
        %final_ind = length(data.lml);

        % output the max lml value
        [lml_max, ind_max] = max(data.lml);
        fprintf("max lml is %d at length=%d and process=%d \n",lml_max,data.length_scale(ind_max),data.process_noise(ind_max));



        figure()
        scatter3(data.length_scale(1:final_ind),data.process_noise(1:final_ind),data.lml(1:final_ind),20,'b.'), hold on
        scatter3(data.length_scale(ind_max),data.process_noise(ind_max),lml_max,'ro','filled')
        xlabel('lengthscale'),ylabel('process noise \sigma_f'),zlabel('LML')
        title(strrep(filename,'_',' '))
    end
end


%%





% load the csv
filename = "systematic_10_lml";
data = readtable(strcat(filename,".csv"));
%data = readtable("exact_100_lml.csv");
% remove rows that have NAN LML
to_remove = isnan(data.lml)';
data = data(~to_remove,:);
data = sortrows(data,'lml','descend');


final_ind = round(0.9*length(data.lml));

% output the max lml value
[lml_max, ind_max] = max(data.lml);
fprintf("max lml is %d at length=%d and process=%d \n",lml_max,data.length_scale(ind_max),data.process_noise(ind_max));


% make a scatter plot
%close all
figure()
scatter3(data.length_scale(1:final_ind),data.process_noise(1:final_ind),data.lml(1:final_ind),10,'b.'), hold on
scatter3(data.length_scale(ind_max),data.process_noise(ind_max),lml_max,'ro','filled')
xlabel('lengthscale'),ylabel('process noise \sigma_f'),zlabel('LML')
title(filename)
