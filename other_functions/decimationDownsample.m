function [X_d, Y_d] = decimationDownsample(X,Y,b)
% INPUTS: X and Y data, and num to downsample by
% If num = 2, will take every second value of X and Y 
% RETURNS: downsampled X and Y as X_d and Y_d

X_d = X(1:b:end);
Y_d = Y(1:b:end);


disp(strcat("Training Data: ", num2str(length(X))," pts."))
disp(strcat("Slim Data: ",num2str(length(X_d))," pts (every 1 in ", num2str(b), "pts kept)."))

end