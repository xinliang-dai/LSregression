function x = line_search_cg(A,b,tol,maxit)
% lin
% check alg options
if nargin < 3 || isempty(tol), tol = 1e-6; end
if nargin < 4 || isempty(maxit), maxit = min(20,length(b)); end
% intial
x = zeros(size(b));
r = A*x -b;
p = - r;
Ap = A*p;
rr = r'*r;
i = 1;
while max(abs(r))>eps &&  i<=maxit
    
    rho = rr/(p'*Ap);       % one-dim minimizer
    x   = x + rho*p;        % update state 
    r   = r + rho*Ap;       % update residual
    rr_new = r'*r;          
    b   = rr_new/rr;        % update the parameter to ensure conjugate 
    rr  = rr_new;
    p   = -r + b*p;
    Ap  = A*p;
    i   = i+1;
end
end