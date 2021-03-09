function x = pcg_steihaug(A,b,tol,maxit,delta,C_T)
% lin
% check alg options
if nargin < 3 || isempty(tol), tol = 1e-6; end
if nargin < 4 || isempty(maxit), maxit = min(20,length(b)); end
if nargin <5 || isempty(delta), delta = 1; end
% intial
tol= tol^2;
x = zeros(size(b));
r = -b;
y = C_T'\(C_T\r);
p = - y;
Ap = A*p;
ry = r'*y;
i = 1;
while max(abs(r))>eps &&  i<=maxit
    if p'*Ap <= 0
        % 1. terminate when negative defined A, return cauchy-point
        tau = find_steplength_on_edge(x,b,delta);
        x   = x + tau * p;
        return 
    end
    rho = ry/(p'*Ap);       % one-dim minimizer
    xk  = x;
    x   = x + rho*p;        % update state 
    if x'*x > delta^2
        % 2. terminate when new step encounters edge of trust-region
        % tau = arg_tau ||x + tau*p||_M = delta
        tau = find_steplength_on_edge(xk,p,delta);
        x   = xk + tau* p;
        return
    end
    r   = r + rho*Ap;       % update residual
    y   = C_T'\(C_T\r);
    ry_new = r'*y;
    if ry_new<tol
        % 3. terminate when reach CG solution within trust region
        return
    end
    beta   = ry_new/ry;        % update the parameter to ensure conjugate 
    ry  = ry_new;
    p   = -y + beta*p;
    Ap  = A*p;
    i   = i+1;
end
end

function tau = find_steplength_on_edge(x,p,delta)
    xx  = x'*x;
    pp  = p'*p;
    xp  = x'*p;
    tau = (-xp + sqrt(xp^2+pp*(delta^2-xx)))/pp; %7.5.5 Trust Region Methods, SIAM
end