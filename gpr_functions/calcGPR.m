function [GPR] = calcGPR(X,Y,hp,data_stats)
% INPUT1,2: training data X, Y 
% INPUT3: original datapts N (not downsampled M)
% OUTPUT1: prediction data x_values X_Star
% OUTPUT2: prediction data y_values Y_Star_Hat
% OUTPUT3: prediction data variances Sigma_Star

% Kernel Hyperparameters as argument
tic

% assign data_stats values to local variable names
nnum = data_stats.nnum;
X_beg = data_stats.beg;
X_end = data_stats.end;

% Prediction Points STAR 
%X_Star = [[(-15+X_beg):2:(15+X_end)]'; X];            % vertical array,x-points we will predict at, INCLUDES INPUT DATA
X_Star = [(-15+X_beg):2:(15+X_end)]';                   % vertical array, x-points we will predict at, EXCLUDES INPUT DATA
                                                        % better for method comparison
 
npts = length(X);                                       % actual pts in our training data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MATRIX CALCS

% Calculate V and Inv(V)                             % depends on training x-points only
W = (hp.sigma_p.^2)*eye(npts);                       % Whitenoise (identity * sigmasquared)
V = K_Function(X,X,hp) + W;                          % Calculate Covariance Matrix using Kernel

% Generate K Parameters
K_Star = K_Function(X_Star,X,hp);                    % Calculate K_Star for New Point(s)
K_StarStar = K_Function(X_Star,X_Star,hp);           % Calculate K_StarStar for New Point(s)

% Cholesky Decomposition
L = chol(V,'lower');

% Calculate Predictions!                             % Finally bring in the training y-points here
Y_Star_Hat = K_Star*CholeskySolve(L,Y);              % Y Predictions (mean values of Gaussians)
CapSigma_Star = K_StarStar-K_Star*CholeskySolve(L,Y);% Variance Predictions (gives us mean, var for each pt)
Y_Star_Var = diag(CapSigma_Star);                    % The diagonals store the variances we want!

% Log Marginal Likelihood
LML = -0.5*log(det(V)) - 0.5*Y'*CholeskySolve(L,Y) - 0.5*nnum*log(2*pi);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUTS

% Update GPR Object with Relevant Output Data
GPR.compute_time = toc;
GPR.X = X;
GPR.Y = Y;
GPR.X_Star = X_Star;
GPR.Y_Star_Hat = Y_Star_Hat;
GPR.Y_Star_Var = Y_Star_Var;
GPR.LML = LML;
GPR.npredict = length(GPR.Y_Star_Hat);
GPR.ntrain = length(GPR.X);


%Plot Training and Prediction Data w/ Error Bars
% figure
% errorbar(X_Star,Y_Star_Hat,Y_Star_Var,'ro'), hold on
% plot(X,Y,'bo','MarkerFaceColor','b'), grid on, hold off
% xlabel('X Values'), ylabel('Y Values')
% legend('Prediction \mu and \sigma^2','Raw Data')

% Bounded Plot (Training Data + Predictions + 2sigma Upper and Lower Bound
sortobj = [X_Star, Y_Star_Hat, Y_Star_Var];
sortobj = sortrows(sortobj);
sorted.X_Star = sortobj(:,1); sorted.Y_Star_Hat = sortobj(:,2); sorted.Y_Star_Var = sortobj(:,3);

figure
p1 = plot(sorted.X_Star,sorted.Y_Star_Hat + 2*sqrt(sorted.Y_Star_Var),'r','LineWidth',2); hold on %upper bound
plot(sorted.X_Star,sorted.Y_Star_Hat - 2*sqrt(sorted.Y_Star_Var),'r','LineWidth',2); % lower bound
p2 = plot(sorted.X_Star,sorted.Y_Star_Hat,'r--','Linewidth',2); % prediction means
p3 = plot(X,Y,'bo','MarkerFaceColor','b','MarkerSize',4); % training data
xlabel('X Values'), ylabel('Y Values'), title('Gaussian Process Regression')
grid on, legend([p3 p2 p1],"Training Data","Prediction \mu","Prediction 2\sigma")

% Update GPR Object with Downsample Info
if npts < nnum % we downsampled!
    GPR.title_str = strcat("GPR with M = ",num2str(npts), " pts");
    GPR.downsampled = 1;
    title(GPR.title_str)
    fprintf("...downsampled GPR (%0.0f pts) compute time: %0.3f\n",npts,GPR.compute_time);
else % we did not downsample
    GPR.title_str = strcat("GPR with N = ",num2str(nnum), " pts");
    GPR.downsampled = 0;
    title(GPR.title_str)
    fprintf("...exact GPR (%0.0f pts) compute time: %0.3f\n",npts,GPR.compute_time);
end


end