function [k] = SqExpKernel(x,y,hp)
% Squared Exponential Kernel Function
% Inputs: two datapoints (x,y) and hyperparameters (hp)

% hyperparameters stored in structure
    % a large sigma_p is a tigther fit - better at/between training data but
    %   blows up beyond the training data


%%%%% EXACT SQUARED EXPONENTIAL KERNEL
% Distance between Points and Covariance Value
if size(x(1,:),2) == 1                              % data is 1x1 scalars
    dist = abs(x-y);                                % 1D
    k = (hp.sigma_p^2)*exp( -(dist^2)/(2*hp.L^2) ); % scalar

elseif size(x(1,:),2) == 2                          % data is 1x2 vectors
    dist = sqrt( (x(1)-y(1))^2 + (x(2)-y(2))^2 );   % 2D
    k = (hp.sigma_p^2)*exp( -(dist^2)/(2*hp.L^2) ); % scalar

else
    disp("Kernel function is NOT setup to process your type of data.")
    return

end






% Great Resource
% http://evelinag.com/Ariadne/covarianceFunctions.html