function [x, flag, relres, iter] = pcgsh(A,b,delta,tol,maxit,M)

% Check the options.
if nargin < 3, delta = nan; end
if nargin < 4 || isempty(tol), tol = 1e-6; end
if nargin < 5 || isempty(maxit), maxit = min(20,length(b)); end
PC = nargin > 5 && (isa(M,'function_handle') || ...
     (isnumeric(M) && all(size(M) == length(b))));

% Initialize PCG-Steihaug.
x = zeros(size(b)); % improve speed
r = -b;             % 
if PC
    if isnumeric(M), y = M\r;
    else y = M(r); end
    d = -y;
    rr = r'*y;
else
    d = -r;
    rr = r'*r;
end
normb = sqrt(rr);  % norm-2
flag = 1;

% PCG-Steihaug.
for iter = 1:maxit

    if isnumeric(A), Ad = A*d;
    else Ad = A(d); end
    alpha = rr/(d'*Ad);

    x1 = x;
    x = x+alpha*d;
    
    % Steihaug's stopping criterion. The case of directions of negative
    % curvature does not need to be handled for NLS.
    if ~isnan(delta) && norm(x) >= delta
        xx = x1'*x1;
        dd = d'*d;
        c = real(x1'*d);
        alpha = (delta^2-xx)/(c+sqrt(c^2+dd*(delta^2-xx)));
        x = x1+alpha*d;
        flag = 2;
    end
    
    r = r+alpha*Ad;
    rr1 = rr;
    if PC
        if isnumeric(M), y = M\r;
        else y = M(r); end
        rr = r'*y;
    else
        rr = r'*r;
    end
    
    if PC, relres = norm(r)/normb;
    else relres = sqrt(rr)/normb; end
    if flag ~= 1, break; end
    if relres < tol, flag = 0; break; end

    beta = rr/rr1;
    if PC, d = -y+beta*d;
    else d = -r+beta*d; end

end

end
