function [X_d, Y_d] = randomDownsample(X,Y,slim_pts)
% Uniform random downsampling of Raw Data to create Slim Data

% Generate unique random indices ranging from 1 to nnum
nnum = length(X);
idcs = zeros(slim_pts,1);

for i = 1:length(idcs)
    new_index = randi(nnum); % uniform random index between 1 and nnum
    while sum(new_index == idcs) > 0 % keep sampling until a unique index found
        new_index = randi(nnum);
    end
    idcs(i) = new_index; % add new index to index list
end

% Generate Slim Data
X_d = X(idcs); Y_d = Y(idcs);

% Output Info
disp(strcat("Training Data: ", num2str(nnum), " pts."))
disp(strcat("Slim Data: ", num2str(slim_pts), " pts."))

end

