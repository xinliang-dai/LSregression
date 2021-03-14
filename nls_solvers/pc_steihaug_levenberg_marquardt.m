function [xsol, flag, logg] = pc_steihaug_levenberg_marquardt(problem,preconditioned)
    % Gauss-Newtond with standard conjugate gradient 
    % intial
    if isempty(problem.options)
        problem.options = nlsOption;
    end
    tol      = problem.options.tol;
    toltol   = tol^2;
    iter_max = problem.options.iter_max;
        
    % create logg from interInfo class
    xk       = problem.x0;
    Nx       = numel(xk);
    logg     = iterInfo(iter_max, Nx);
    
    % initial mat/vec for alg
    r        = problem.r;            % residual fucntion
    Jk       = problem.drdx;         % jacobian function
    r_vector = r(xk);                % initial residual vector   
    fval_k   = r_vector'*r_vector/2; % initial cost
    % initial sensitivities of cost
    [grad, B_mat] = update_sensitivities(Jk,xk,r_vector);    
    % preconditioning
    if preconditioned
        dk            = diag(B_mat);
        dk(dk==0)     = tol;
        Dk            = sparse(diag(dk));
    else
        Dk            = speye(Nx);
    end

    
    % select linsolver
    switch problem.options.cg_method
        case {'without CG'}
            cg_steihaug_activated=false;
            delta    = [];
        case {'CG-Steihaug'}
            cg_steihaug_activated=true;
            % initial trust region for cg-steihaug
            eta2     = 0.9;                % quadratic model very accurate?
            theta1   = 2.5;                % multiplier for increasing radius of trust region
            theta2   = 0.25;               % multiplier for decreasing radius of trust region
            delta    = sqrt(grad'*grad);     % initial radius of trust region            
        otherwise 
            warning('CG method not found, using CG-Steihaug method')
            cg_steihaug_activated=true;
            % initial trust region for cg-steihaug
            eta2     = 0.9;                % quadratic model very accurate?
            theta1   = 2.5;                % multiplier for increasing radius of trust region
            theta2   = 0.25;               % multiplier for decreasing radius of trust region
            delta    = sqrt(grad'*grad);     % initial radius of trust region
    end
    
    % Nielsen damping strategies of Levenberg-Marquardt method
    v        = 2;
    mu       = max(abs(diag(B_mat)));
    
    % initial success_flag and while loop
    flag     = false;
    i        = 1;
    
    % start Levenberg iterates
    while i<=iter_max && ~flag
        % Levenberg-Marquardt step
        if cg_steihaug_activated
            pk        = cg_steihaug(B_mat+mu*Dk,-grad,toltol,Nx,delta);
        else
            pk        = - [B_mat+mu*Dk]\grad;
        end
        pp            = sqrt(pk'*pk);
        r_vector      = r(xk+pk);                         % residual vector   
        logg.fval(i)  = r_vector'*r_vector/2;     % cost value    
        logg.dfval(i) = fval_k-logg.fval(i);              % decreasing of cost   
        % terminate if steplength is small enough
        if  pp<=tol*(sqrt(xk'*xk)+tol) % || abs(logg.dfval(i))<=tol
            flag = true;
        end
        % updating damping parameter
        dmval         = 0.5*pk'*(mu*Dk*pk-grad);  % mk(xk)-mk(xk+pk)  
        rho           = logg.dfval(i)/dmval;          % trustworthness rho
        % updating trustworthness rho and radius of trust-region
        if rho>0
            % descent direction, pk accepted
            xk = xk+pk;
            % decreasing lambda -> approaching GN method
            mu        = mu * max(1/3, 1-(2*rho-1)^3);
            v         = 2;
            fval_k    = logg.fval(i);       % record new cost            
            % updating sensitivities
            [grad, B_mat] = update_sensitivities(Jk,xk,r_vector);
            % preconditioning
            if preconditioned
                dk            = diag(B_mat);
                dk(dk==0)     = tol;
                Dk            = diag(dk);
            end
            % terminate if first order derivative small enought
            if  abs(max(grad))<=tol
                flag  = true;
            end
            % increase trust region for CG-Steihaug method, when good fit
            if   cg_steihaug_activated && rho>eta2
                delta = max(theta1*sqrt(pp),delta);
            end
        else
            % not descont direction, reject pk
            % increasing lambda -> approaching steepest descent
            mu    = mu*2;
            v     = v * 2;
            pk    = 0;            
            % reduce radius of trust region for CG-Steihaug method, when
            % bad fit
            if cg_steihaug_activated
                delta = theta2*sqrt(pp);
            end
        end

        % recording iter Info of current iteration
        if cg_steihaug_activated
            logg.delta(i) = delta;        
        end
        if cg_steihaug_activated
            logg.delta(i) = delta;        
        end        
        logg.dmval(i) = dmval;
        logg.rho(i)   = rho;
        logg.mu(i)    = mu;
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

function x = cg_steihaug(A,b,toltol,maxit,delta)
% lin
% check alg options
if nargin < 3 || isempty(toltol), toltol = 1e-12; end
if nargin < 4 || isempty(maxit), maxit = min(20,length(b)); end
% intial
x = zeros(size(b));
r = -b;  
% without preconditioning
p = -r;
rr = r'*r;
Ap = A*p;
i = 0;
while max(abs(r))>eps &&  i<=maxit
    % A always positive defined
    if ~isempty(delta) && p'*Ap <= 0  
        % 1. terminate when negative defined A, return cauchy-point
        tau = find_steplength_on_edge(x,p,delta);
        x   = x + tau * p;
        return 
    end
    rho = rr/(p'*Ap);       % one-dim minimizer
    xk  = x;
    x   = x + rho*p;        % update state
    if ~isempty(delta) && x'*x > delta^2 
        % 2. terminate when new step encounters edge of trust-region
        tau = find_steplength_on_edge(xk,p,delta);
        x   = xk + tau * p;
        return
    end
    r      = r + rho*Ap;       % update residual
    rr_new = r'*r;
    if rr_new<toltol
        % 3. terminate when reach CG solution within trust region
        return
    end
    beta   = rr_new/rr;      % update the parameter to ensure conjugate 
    rr  = rr_new;
    p   = -r + beta*p;
    Ap  = A*p;
    i   = i+1;
end
end

function tau = find_steplength_on_edge(x,p,delta)
    % tau = arg_tau ||x + tau*p||_M = delta
    xx  = x'*x;
    pp  = p'*p;
    xp  = x'*p;
    tau = (-xp + sqrt(xp^2+pp*(delta^2-xx)))/pp; %7.5.5 Trust Region Methods, SIAM
end