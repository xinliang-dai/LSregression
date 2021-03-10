function x = cg_steihaug(A,b,tol,maxit,delta)
% lin
% check alg options
if nargin < 3 || isempty(tol), tol = 1e-6; end
if nargin < 4 || isempty(maxit), maxit = min(20,length(b)); end
if nargin <5 || isempty(delta), delta = 1; end
delta=delta^2;
% intial
tol= tol^2;
x = zeros(size(b));
r =-b;
p = - r;
Ap = A*p;
rr = r'*r;
i = 1;
while max(abs(r))>eps &&  i<=maxit
    if p'*Ap <= 0
        % 1. terminate when negative defined A, return cauchy-point
        tau = find_steplength_on_edge(x,b,delta);
        x   = x + tau * p;
        return 
    end
    rho = rr/(p'*Ap);       % one-dim minimizer
    xk  = x;
    x   = x + rho*p;        % update state 
    if x'*x > delta
        % 2. terminate when new step encounters edge of trust-region
        tau = find_steplength_on_edge(xk,p,delta);
        x   = xk + tau * p;
        return
    end
    r   = r + rho*Ap;       % update residual
    rr_new = r'*r;
    if rr_new<tol
        % 3. terminate when reach CG solution within trust region
        return
    end
    b   = rr_new/rr;        % update the parameter to ensure conjugate 
    rr  = rr_new;
    p   = -r + b*p;
    Ap  = A*p;
    i   = i+1;
end
end

function tau = find_steplength_on_edge(x,p,delta)
    xx  = x'*x;
    pp  = p'*p;
    xp  = x'*p;
    tau = (-xp + sqrt(xp^2+pp*(delta-xx)))/pp; %7.5.5 Trust Region Methods, SIAM
end