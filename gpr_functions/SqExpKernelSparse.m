function [k] = SqExpKernelSparse(x,y,hp)
% Squared Exponential Kernel Function
% Inputs: two datapoints (x,y) and hyperparameters (hp)

% hyperparameters stored in structure
    % a large sigma_p is a tigther fit - better at/between training data but
    %   blows up beyond the training data


%%%%% APPROXIMATE SQUARED EXPONENTIAL KERNEL

% Calculate Distance between Points
if size(x(1,:),2) == 1                              % data is 1x1 scalars
    dist = abs(x-y);                                % 1D

elseif size(x(1,:),2) == 2                          % data is 1x2 vectors
    dist = sqrt( (x(1)-y(1))^2 + (x(2)-y(2))^2 );   % 2D

else                                                % other
    disp("Kernel function is NOT setup to process your type of data.")
    return

end

% Calculate Covariance Value
if dist > hp.L                                      % Dist > Lengthscale, then 0
    k = 0;
else                                                % Dist < Lengthscale, kernel value
    k = hp.sigma_p^2 * ( (2+cos(2*pi*dist/hp.L))/3*(1-dist/hp.L) + 1/(2*pi)*sin(2*pi*dist/hp.L)  );
end

% Optionally Round Covariance
%k = round(k,4);





% Great Resource
% http://evelinag.com/Ariadne/covarianceFunctions.html