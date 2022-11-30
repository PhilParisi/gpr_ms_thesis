function [LML] = calcLML(L,Y,nnum)
% Calculate the Log Marginal Likelihood (LML)
    % input, CholeskyFactor L, Training Data Y's, nnum total pts N
    % this method uses the cholesky factor

% Note: MATLAB's log function is a natural log (which matches GPR LML equations)
% Without Cholesky: LML = -0.5*log(det(V)) - 0.5*Y'*CholeskySolve(L,Y) - 0.5*nnum*log(2*pi);

% Log Marginal Likelihood Forumla (w/ Cholesky)
LML = -log(prod(diag(L),'all')) - 0.5*Y'*CholeskySolve(L,Y) - 0.5*nnum*log(2*pi);

end

