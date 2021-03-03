function [xsol, flag, logg] = basic_gauss_newton(problem)
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
    df = problem.dfdx(x0);
    f  = problem.f(x0);
    flag  = false;
    i =1;
    while i<iter_max && ~flag
        % gauss-newton step
        p  = - (df'*f)/(df'*df);
        % relative steplength
        rel_steplength = norm(p)/norm(x0);
        if isnan(rel_steplength), rel_steplength = 0; end
        df = problem.dfdx(x0+p);
        f  = problem.f(x0+p);
        df_val = norm(df,inf);
        f_val  = norm(f,inf);
%         norm(dfval)
        x0 = x0+p;
        if rel_steplength<=tol
            xsol = x0;
            flag = true;
        end
        % recording in current iteration
        logg.xk(:,i)           = x0;
        logg.iter(i)           = i;
        logg.fval(i)           = f_val;
        logg.dfval(i)          = df_val;
        logg.rel_steplength(i) = rel_steplength;
        i = i + 1;
    end
    if isnan(xsol)
        xsol = x0;
    end
end