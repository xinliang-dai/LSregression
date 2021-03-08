function x = preconditioned_cg(A,b,tol,maxit,C_T)
% check alg options
if nargin < 3 || isempty(tol), tol = 1e-6; end
norm_gradient = norm(b,2);
tol =min(0.5,norm_gradient)*norm_gradient;
if nargin < 4 || isempty(maxit), maxit = min(20,length(b)); end
% intial
x = zeros(size(b));
r = A*x -b;
y = C_T'\(C_T\b);
p = - y;
Ap = A*p;
ry = r'*y;
i = 1;

while sum(abs(r))>tol &&  i<=maxit
rho = ry/(p'*Ap);       % one-dim minimizer
x   = x + rho*p;        % update state 
r   = r + rho*Ap;       % update residual
y   = C_T'\(C_T\r);    
ry_new = r'*y;
b   = ry_new/ry;        % update the parameter to ensure conjugate 
ry  = ry_new;
p   = -y + b*p;
Ap  = A*p;
i   = i+1;
end
end