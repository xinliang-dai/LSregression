function [xsol, flag, logg] = doglet_gauss_newton(problem)
    % Gauss-Newtond with standard conjugate gradient 
    % intial
    if isempty(problem.options)
        problem.options = nlsProblem;
    else
        tol      = problem.options.tol;
        iter_max = problem.options.iter_max;
    end
    toltol      = tol^2;    
    
    % setting for select CG-method
    if ismember(problem.options.cg_method,'without CG')
        % using inverse to compute Gauss-Newton step
        precondition = false;
        CG_activated = false;
    else
        CG_activated = true;
        switch problem.options.cg_method
            case {'standard CG'}
                precondition = false;
            case {'preconditioned CG'}
                precondition = true;
%             case {'CG-Steihaug'}
%                 trust_region = true;
%                 precondition = false;
%             case {'preconditioned CG-Steihaug'} 
%                 trust_region = true;
%                 precondition = true;
            otherwise
                warning('CG method not found, using CG method')
                precondition = true;
        end
    end
    
    % initial Gauss-Newton method
    xk       = problem.x0;
    xx0      = xk'*xk;
    r        = problem.r;            % residual fucntion
    Jk       = problem.drdx;         % jacobian function
    r_vector = r(xk);                % initial residual vector   
    Jk_mat   = sparse(Jk(xk));       % initial jacobian matrix
    B_mat    = Jk_mat'*Jk_mat;       % initial hessian approximation for LS
    grad     = Jk_mat'*r_vector;     % initial gradient
    fval_k   = r_vector'*r_vector/2; % initial cost
    
    % create logg from interInfo class
    Nx       = numel(xk);
    logg     = iterInfo(iter_max, Nx);
    
    % initial success_flag and while loop
    flag     = false;
    i        = 1;
    
    % initial preconditioning
    if precondition
        options.diagcomp = 10^-6;
    end
    
    % initial trust region - typical choices described by "Trust Region Methods"
    eta1     = 0.05;               % trial step acceptable?
    eta2     = 0.9;                % quadratic model very accurate?
    theta1   = 2.5;                % multiplier for increasing radius of trust region
    theta2   = 0.25;               % multiplier for decreasing radius of trust region
    delta    = sqrt(grad'*grad);     % initial radius of trust region
    
    % start Gauss-Newton iterates
    while i<=iter_max && ~flag
        dd       = delta^2;
        if precondition
            C_T      = ichol(B_mat,options);        % incomplete cholosky to determine M = C^T * C;
        else
            C_T      = [];
        end        
        % Gauss-Newton step
        if ~CG_activated
            % Gauss-Newton without CG method
            p_gn = - B_mat\ grad;
        else
            % Gauss-Newton with CG method
            p_gn = preconditioned_cg(B_mat,-grad,toltol,iter_max,C_T);
        end
        pp_gn         = p_gn'*p_gn;
        if pp_gn <=dd
            pk = p_gn;
        else
            p_sd = - ((grad'*grad) / (grad'*B_mat*grad))* grad; % steepest descent direction
            pp_sd = p_sd'*p_sd;
            if pp_sd >= dd
                pk = sqrt(dd/pp_sd)*p_sd;
            else
                dp = p_gn-p_sd;
                c  = p_sd'*dp;
                if c<=0
                    beta = (-c + sqrt(c^2+dp'*dp*(dd-pp_sd))) / dp'*dp;
                else
                    beta = (dd-pp_sd)/(c+sqrt(c^2+dp'*dp*(dd-pp_sd)));
                end
                pk  = p_sd + beta * dp;
            end
        end
        pp            = pk'*pk;
        r_vector      = r(xk+pk);                         % residual vector   
        logg.fval(i)  = r_vector'*r_vector/2;         
        logg.dfval(i) = fval_k-logg.fval(i);              % cost-value    
        % terminate if steplength is small enough
        if  pp<=tol*(sqrt(xk'*xk)+tol) % || abs(logg.dfval(i))<=tol
            flag = true;
        end
        % updating radius of trust region
        
        dmval         = - pk'*grad - pk'*B_mat*pk/2;  % mk(xk)-mk(xk+pk)  
        rho           = logg.dfval(i)/dmval;          % trustworthness rho
        % updating trustworthness rho and radius of trust-region
        if rho<=eta1  
            % trial step not acceptable, shrink size of TR
            delta         = theta2*sqrt(pp);
            pk            = 0;
        else
            % trial step accaptable
            xk = xk+pk;
            % approxiamation is very satisfying, expand size of TR
            if rho>eta2
                delta     = max(theta1*sqrt(pp),delta);
            end
            % updating sensitivities
            [grad, B_mat] = update_sensitivities(Jk,xk,r_vector);
            fval_k        = logg.fval(i);       % record new cost
        end

        % recording iter Info of current iteration
        logg.delta(i) = delta; 
        logg.dmval(i) = dmval;
        logg.rho(i)   = rho;
        logg.pk(:,i)  = pk;
        logg.grad(:,i)= grad;
        logg.xk(:,i)  = xk;
        logg.iter(i)  = i;
        i = i + 1;
    end
    xsol = xk;
end

function [grad,B_mat] = update_sensitivities(Jk,xk,r_vector)
    Jk_mat  = sparse(Jk(xk));
    grad    = Jk_mat'*r_vector;
    B_mat   = Jk_mat'*Jk_mat;
end

function x = preconditioned_cg(A,b,toltol,maxit,C_T)
% lin
% check alg options
if nargin < 3 || isempty(toltol), toltol = 1e-12; end
if nargin < 4 || isempty(maxit), maxit = min(20,length(b)); end
% if nargin <5 || isempty(delta), delta = 1; end
% intial
x = zeros(size(b));
r = -b;
if ~isempty(C_T)
    % preconditioning activated
    y = C_T'\(C_T\r);
    p = - y;
    rr = r'*y;    
else
    % without preconditioning
    p = -r;
    rr = r'*r;
end
Ap = A*p;
i = 1;
while max(abs(r))>eps &&  i<=maxit
    rho = rr/(p'*Ap);       % one-dim minimizer
    x   = x + rho*p;        % update state
    r   = r + rho*Ap;       % update residual
    if ~isempty(C_T)
        % preconditioning activated
        y   = C_T'\(C_T\r);
        rr_new = r'*y;
    else
        % without preconditioning
        rr_new = r'*r;
    end
    if rr_new<toltol
        % 3. terminate when reach CG solution within trust region
        return
    end
    beta   = rr_new/rr;      % update the parameter to ensure conjugate 
    rr  = rr_new;
    if ~isempty(C_T)
        % preconditioning activated
        p   = -y + beta*p;
    else
        % without preconditioning
        p   = -r + beta*p;
    end    
    Ap  = A*p;
    i   = i+1;
end
end

function tau = find_steplength_on_edge(x,p,dd)
    % tau = arg_tau ||x + tau*p||_M = delta
    xx  = x'*x;
    pp  = p'*p;
    xp  = x'*p;
    tau = (-xp + sqrt(xp^2+pp*(dd-xx)))/pp; %7.5.5 Trust Region Methods, SIAM
end