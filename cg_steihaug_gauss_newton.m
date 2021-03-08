function [xsol, flag, logg] = cg_steihaug_gauss_newton(problem)
    % Gauss-Newtond with standard conjugate gradient 
    % intial
    if isempty(problem.options)
        problem.options = nlsProblem;
    else
        tol      = problem.options.tol;
        iter_max = problem.options.iter_max;
    end
    x0       = problem.x0;
    Nx       = numel(x0);

    % create logg from interInfo class
    logg = iterInfo(iter_max, Nx); 
    
    % start Gauss-Newton method
    r      = problem.r;          % residual fucntion
    Jk     = problem.drdx;       % jacobian function
    r_val  = r(x0);              % initial residual
    Jk_val = Jk(x0);             % initial jacobian matrix
    Jk_norm = norm(Jk_val,2);
    B      = @(x) (Jk(x))'*Jk(x);% approximation of Hessian
    flag  = false;
    i =1;
    while i<iter_max && ~flag 
        B_mat = B(x0);
        % gauss-newton step
        p_class  = - B_mat\ (Jk_val'*r_val);
        % gauss-newton step by conjugate gradient
        p     = cg_steihaug(B_mat,-Jk_val'*r_val,tol,iter_max,delta);
        % updating radius of trust-region
        
        % relative steplength
        rel_steplength = norm(p,2)/norm(x0,2);
        if isnan(rel_steplength), rel_steplength = 0; end
        Jk_val     = Jk(x0+p);
        r_val      = problem.r(x0+p);
        cost       = norm(r_val,2)/2;        
        Jk_norm_new   = norm(Jk_val, 2);
        rho        = (Jk_norm_new - Jk_norm)/Jk_norm;
        Jk_norm    = Jk_norm_new;
%         norm(dfval)
        x0 = x0+p;
        
        if  rel_steplength<=tol 
            xsol = x0;
            flag = true;
        end
        % recording in current iteration
        logg.xk(:,i)           = x0;
        logg.iter(i)           = i;
        logg.fval(i)           = cost;
        logg.dfval(i)          = Jk_norm;
        logg.rel_steplength(i) = rel_steplength;
        i = i + 1;
    end
        xsol = x0;
end