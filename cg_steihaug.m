function x = cg_steihaug(A,b,tol,maxit,delta)
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
    if p'*Ap <= 0
        % negative defined A, return cauchy-point
        tau = find_steplength_on_tr(x,p,delta);
        x   = x + tau * p;
        return 
    end
    rho = rr/(p'*Ap);       % one-dim minimizer
    x0  = x;
    x   = x + rho*p;        % update state 
    if norm(x,2)>delta
        % encounter edge of trust-region
        tau = find_steplength_on_tr(x0,p,delta);
        x   = x0 + tau * p;
        return
    end
    r   = r + rho*Ap;       % update residual
    rr_new = r'*r;          
    b   = rr_new/rr;        % update the parameter to ensure conjugate 
    rr  = rr_new;
    p   = -r + b*p;
    Ap  = A*p;
    i   = i+1;
end
end

function tau = find_steplength_on_tr(x,p,delta)
  %  tau = (delta^2-x'*x-p'*p)/ (2*p'*x);
end