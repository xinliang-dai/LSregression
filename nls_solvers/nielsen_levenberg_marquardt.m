function [xsol, flag, logg] = nielsen_levenberg_marquardt(problem)
    % Gauss-Newtond with standard conjugate gradient 
    % intial
    if isempty(problem.options)
        problem.options = nlsProblem;
    else
        tol      = problem.options.tol;
        iter_max = problem.options.iter_max;
    end

    
    xk       = problem.x0;
    xx0      = xk'*xk;

    toltol   = tol^2;
    
    % create logg from interInfo class
    Nx       = numel(xk);
    logg     = iterInfo(iter_max, Nx);
    
    r        = problem.r;            % residual fucntion
    Jk       = problem.drdx;         % jacobian function
    r_vector = r(xk);                % initial residual vector   
    Jk_mat   = sparse(Jk(xk));       % initial jacobian matrix
    B_mat    = Jk_mat'*Jk_mat;       % initial hessian approximation for LS
    grad     = Jk_mat'*r_vector;     % initial gradient
    fval_k   = r_vector'*r_vector/2; % initial cost
    dk       = zeros(Nx,1);
    Dk       = zeros(Nx,Nx);
    % initial Levenberg-Marquadt method
%     eta1     = 0.05;               % trial step acceptable?
%     eta2     = 0.9;                % quadratic model very accurate?
%     theta1   = 2.5;                % multiplier for increasing radius of trust region
%     theta2   = 0.25;               % multiplier for decreasing radius of trust region
%     lambda   = max(abs(grad));     % initial radius of trust region
    eye_null = speye(Nx,Nx);
    vec_null = sparse(zeros(Nx,1));    
    
    % Nielsen damping strategies of Levenberg-Marquardt method
    v        = 2;
    lambda   = max(diag(B_mat));
    
    % initial success_flag and while loop
    flag     = false;
    i        = 1;
    
    % start Gauss-Newton iterates
    while i<=iter_max && ~flag
        % Levenberg-Marquardt step
        dk            = diag(B_mat);
        dk(dk==0)      = tol;
        Dk            = diag(dk);        
        pk            = -[B_mat+lambda*Dk]\grad;
        pp            = sqrt(pk'*pk);        
        r_vector_new  = r(xk+pk);                         % residual vector   
        logg.fval(i)  = r_vector_new'*r_vector_new/2;         
        logg.dfval(i) = fval_k-logg.fval(i);              % cost-value    
        % terminate if steplength is small enough
        if  max(abs(pk))<=tol*(max(xk)+tol)  || abs(logg.dfval(i))<=tol
            flag = true;
        end
        % updating damping parameter
        dmval         = 0.5*pk'*(lambda*Dk*pk-grad);  % mk(xk)-mk(xk+pk)  
        rho           = logg.dfval(i)/dmval;          % trustworthness rho
        % updating trustworthness rho and radius of trust-region
        if rho>0
            % descent direction, pk accepted
            xk = xk+pk;
            % decreasing lambda -> approaching GN method
            lambda        = lambda * max(1/3, 1-(2*rho-1)^3);
            v             = 2;
            % updating sensitivities
            [grad, B_mat] = update_sensitivities(Jk,xk,r_vector_new);
            if  abs(max(grad))<=tol
                flag = true;
            end           
            fval_k        = logg.fval(i);       % record new cost            
        else
            % not descont direction, reject pk
            % increasing lambda -> approaching steepest descent
            mu = mu*2;
            v  = v * 2;
        end

        % recording iter Info of current iteration
        logg.dmval(i) = dmval;
        logg.rho(i)   = rho;
        logg.delta(i) = lambda; 
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