function x = conjugate_gradient(A,b,tol,maxit)
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
while norm(r,2)>tol &&  i<=maxit
    rho = rr/(p'*Ap);
    x   = x + rho*p;
    r   = r + rho*Ap;
    rr_new = r'*r;
    b   = rr_new/rr;
    rr  = rr_new;
    p   = -r + b*p;
    i   = i+1;   
end
end