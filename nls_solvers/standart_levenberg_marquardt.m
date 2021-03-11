function [xsol, flag, logg] = standart_levenberg_marquardt(problem)
    % Gauss-Newtond with standard conjugate gradient 
    % intial
    if isempty(problem.options)
        problem.options = nlsProblem;
    else
        tol      = problem.options.tol;
        iter_max = problem.options.iter_max;
    end
    toltol      = tol^2;    
    
%     % setting for select CG-method
%     if ismember(problem.options.cg_method,'without CG')
%         % using inverse to compute Gauss-Newton step
%         CG_activated = false;  
%         trust_region = false;
%         precondition = false;
%     else
%         CG_activated = true;
%         switch problem.options.cg_method
%             case {'standard CG'}
%                 trust_region = false;
%                 precondition = false;
%             case {'preconditioned CG'}
%                 trust_region = false;
%                 precondition = true;
%             case {'CG-Steihaug'}
%                 trust_region = true;
%                 precondition = false;
%             case {'preconditioned CG-Steihaug'} 
%                 trust_region = true;
%                 precondition = true;
%             otherwise
%                 warning('CG method not found, using preconditioned CG-Steihaug method')
%                 trust_region = true;
%                 precondition = true;
%         end
%     end
    
    xk       = problem.x0;
    xx0      = xk'*xk;
    
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

    % initial Levenberg-Marquadt method
    eta1     = 0.05;               % trial step acceptable?
    eta2     = 0.9;                % quadratic model very accurate?
    theta1   = 2.5;                % multiplier for increasing radius of trust region
    theta2   = 0.25;               % multiplier for decreasing radius of trust region
    lambda   = max(abs(grad));     % initial radius of trust region
    eye_null = speye(Nx,Nx);
    vec_null = sparse(zeros(Nx,1));    
    

    
    % initial success_flag and while loop
    flag     = false;
    i        = 1;
    
%     % initial preconditioning
%     if precondition
%         C_T      = ichol(B_mat);        % incomplete cholosky to determine M = C^T * C;
%     else
%         C_T      = [];
%     end
%     
%     % initial trust region - typical choices described by "Trust Region Methods"
%     if trust_region

%     else
%         delta    = [];
%     end
    
    % start Gauss-Newton iterates
    while i<=iter_max && ~flag
%         B_lm = [B_mat, sqrt()+]
        D          = ichol(B_mat);
        % Levenberg-Marquardt step
        pk            = -[B_mat+lambda*eye_null]\grad;
        pp            = pk'*pk;        
        r_vector_new  = r(xk+pk);                         % residual vector   
        logg.fval(i)  = r_vector_new'*r_vector_new/2;         
        logg.dfval(i) = fval_k-logg.fval(i);              % cost-value    
        % terminate if steplength is small enough
        if  pp<=toltol*xx0  || abs(logg.dfval(i))<=tol
            flag = true;
        end
        % updating radius of trust region
%         if trust_region
            dmval         = - pk'*grad - pk'*B_mat*pk/2;  % mk(xk)-mk(xk+pk)  
            rho           = logg.dfval(i)/dmval;          % trustworthness rho
            % updating trustworthness rho and radius of trust-region
            if rho<=eta1  
                % trial step not acceptable, shrink size of TR
                lambda         = theta2*sqrt(pp);
                pk            = 0;
            else
                % trial step accaptable
                xk = xk+pk;
                % approxiamation is very satisfying, expand size of TR
                if rho>eta2
                    lambda     = max(theta1*sqrt(pp),lambda);
                end
                % updating sensitivities
                [grad, B_mat] = update_sensitivities(Jk,xk,r_vector_new);
                fval_k        = logg.fval(i);       % record new cost
            end
            logg.delta(i) = lambda; 
%         else
%             % Trust Region not activated
%             xk            = xk+pk;
%             % updating sensitivities
%             [grad, B_mat] = update_sensitivities(Jk,xk,r_vector_new);
%             fval_k        = logg.fval(i);       % new cost
%         end
        % recording iter Info of current iteration
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
% 
% function x = pcg_steihaug(A,b,toltol,maxit,delta,C_T)
% % lin
% % check alg options
% if nargin < 3 || isempty(toltol), toltol = 1e-12; end
% if nargin < 4 || isempty(maxit), maxit = min(20,length(b)); end
% % if nargin <5 || isempty(delta), delta = 1; end
% % intial
% x = zeros(size(b));
% r = -b;
% if ~isempty(C_T)
%     % preconditioning activated
%     y = C_T'\(C_T\r);
%     p = - y;
%     rr = r'*y;    
% else
%     % without preconditioning
%     p = -r;
%     rr = r'*r;
% end
% Ap = A*p;
% i = 1;
% while max(abs(r))>eps &&  i<=maxit
%     if ~isempty(delta) && p'*Ap <= 0  
%         % 1. terminate when negative defined A, return cauchy-point
%         tau = find_steplength_on_edge(x,b,delta);
%         x   = x + tau * p;
%         return 
%     end
%     rho = rr/(p'*Ap);       % one-dim minimizer
%     xk  = x;
%     x   = x + rho*p;        % update state
%     if ~isempty(delta) && x'*x > delta^2 
%         % 2. terminate when new step encounters edge of trust-region
%         tau = find_steplength_on_edge(xk,p,delta);
%         x   = xk + tau * p;
%         return
%     end
%     r   = r + rho*Ap;       % update residual
%     if ~isempty(C_T)
%         % preconditioning activated
%         y   = C_T'\(C_T\r);
%         rr_new = r'*y;
%     else
%         % without preconditioning
%         rr_new = r'*r;
%     end
%     if rr_new<toltol
%         % 3. terminate when reach CG solution within trust region
%         return
%     end
%     beta   = rr_new/rr;      % update the parameter to ensure conjugate 
%     rr  = rr_new;
%     if ~isempty(C_T)
%         % preconditioning activated
%         p   = -y + beta*p;
%     else
%         % without preconditioning
%         p   = -r + beta*p;
%     end    
%     Ap  = A*p;
%     i   = i+1;
% end
% % end
% 
% function tau = find_steplength_on_edge(x,p,delta)
%     % tau = arg_tau ||x + tau*p||_M = delta
%     xx  = x'*x;
%     pp  = p'*p;
%     xp  = x'*p;
%     tau = (-xp + sqrt(xp^2+pp*(delta^2-xx)))/pp; %7.5.5 Trust Region Methods, SIAM
% end