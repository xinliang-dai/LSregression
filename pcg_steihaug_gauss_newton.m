function [xsol, flag, logg] = pcg_steihaug_gauss_newton(problem)
    % Gauss-Newtond with standard conjugate gradient 
    % intial
    if isempty(problem.options)
        problem.options = nlsProblem;
    else
        tol      = problem.options.tol;
        iter_max = problem.options.iter_max;
    end
    tol     = tol^2;    
    
    % initial Gauss-Newton method
    xk       = problem.x0;
    xx0      = xk'*xk;
    r        = problem.r;          % residual fucntion
    Jk       = problem.drdx;       % jacobian function
    
    r_vector = r(xk);              % initial residual vector   
    Jk_mat   = sparse(Jk(xk));             % initial jacobian matrix
    B_mat    = Jk_mat'*Jk_mat;     % initial hessian approximation for LS
    grad     = Jk_mat'*r_vector;   % initial gradient
    fval_k   = r_vector'*r_vector/2; % initial cost
    
    % create logg from interInfo class
    Nx       = numel(xk);
    logg     = iterInfo(iter_max, Nx); 
    
    % initial trust region 
    eta1     = 0.05;               % trial step acceptable?
    eta2     = 0.9;                % quadratic model very accurate?
    theta1   = 2.5;                % multiplier for increasing radius of trust region
    theta2   = 0.25;               % multiplier for decreasing radius of trust region
    delta    = max(abs(grad));      % initial radius of trust region
    
    
    flag     = false;
    i        = 1;
    while i<=iter_max && ~flag 

        
        % gauss-newton step by preconditioned CG-Steihaug
        C_T          = ichol(B_mat);
        p_pcg        = pcg_steihaug(B_mat,-grad,tol,iter_max,delta,C_T);
%         p_p        = preconditioned_cg(B_mat,-grad,tol,iter_max,C_T);
%         p            = cg_steihaug(B_mat,-grad,tol,iter_max,delta);
%          norm(p-p_pcg,2)
        % updating trustworthness of trust-region
        logg.dfval(i) = - p_pcg'*grad - p_pcg'*B_mat*p_pcg/2;   % mk(xk)-mk(xk+pk)      
        r_vector_new  = r(xk+p_pcg);                           % residual vector   
        logg.fval(i)  = r_vector_new'*r_vector_new/2;         % cost-value    
        rho           = (fval_k-logg.fval(i))/logg.dfval(i);  % trustworthness rho
        pp            = p_pcg'*p_pcg;
        % updating radius of trust-region
        if rho<=eta1
            
            % trial step not acceptable, shrink size of TR
            delta         = theta2*sqrt(pp);
            
        else
            
            % trial step accaptable
            xk = xk+p_pcg;
            logg.pk(:,i)  = p_pcg;
            
            % approxiamation is very satisfying, expand size of TR
            if rho>eta2
                delta     = max(theta1*sqrt(pp),delta);
            end        
            
            % updating parameters    
            r_vector      = r_vector_new;       % new residual vector   
            Jk_mat        = sparse(Jk(xk));             % new jacobian matrix
            B_mat         = Jk_mat'*Jk_mat;     % new hessian approximation for LS
            grad          = Jk_mat'*r_vector;   % new gradient
            fval_k        = logg.fval(i);       % new cost
            
            % terminate if termination condition satisfied
            if  pp<=tol*xx0 || grad'*grad<=tol
                xsol = xk;
                flag = true;
                return;
            end
        end
        % recording in current iteration
        logg.xk(:,i)  = xk;
        logg.iter(i)  = i;
        logg.delta(i) = delta;  
        i = i + 1;
    end
        xsol = xk;
end