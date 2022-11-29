function x = CholeskySolve(L,b)
% solve Ax = b using cholesky decomposition
% input is L, cholesky factor of matrix that needs to be inverted
% input is b, matrix to multiply inverse against
% note that the cholesky factor needs to have error added to it!!

% Ax = b --> Lx = b per notation. L is cholesky factor chol(J,'lower')
% essentially, if you have U*V^-1*W
% you should L = chol(V,'lower'); then replace w/ U*CholeskySolve(L,W)

    x = (L')\(L\b);

end