classdef nlsProblem
    %   define standard nonlinear problem
    %   the nonlinear least squares objective function F(x) := 0.5*(f(x)'*f(x)).
    properties
        % residual
        r 
        % jacobian matrix - first derivative of residual
        drdx 
        % initial point
        x0         double    {mustBeNumeric}
        %
        bound      double    {mustBeNumeric}
        options    nlsOption
    end
    
    methods
        % conductor
        function obj = nlsProblem(r,drdx,x0,bound)
            % residual vector-function
            obj.r    = r;
            % jacobian matrix
            obj.drdx = drdx;
            % initial point
            obj.x0   = x0;
            % lower and upper bounds
            if nargin >3 && ~isempty(bound)
                obj.bound = bound;
                obj.check_bound_dim;
            end
        end        
        % solve the nonlinear least-squares problem
        function [xsol, logg, flag, t]= solve_nls(obj)
            % switch from different methods
            if isempty(obj.options)
                obj.options = nlsOption;
            end
            switch obj.options.nls_method
                case {'Gauss-Newton'}
                    tic
                    [xsol,flag,logg] = pcg_steihaug_gauss_newton(obj);
                    t = toc;
                    logg = logg.iter_dataprocessing;
                case {'Levenberg-Marquardt'}
                    preconditioned   = false;
                    tic
                    [xsol,flag,logg] = pc_steihaug_levenberg_marquardt(obj,preconditioned);
                    t = toc;
                    logg             = logg.iter_dataprocessing;
                case {'preconditioned Levenberg-Marquardt'}
                    preconditioned   = true;
                    tic
                    [xsol,flag,logg] = pc_steihaug_levenberg_marquardt(obj,preconditioned);
                    t = toc;
                    logg = logg.iter_dataprocessing;                    
                case {'Dogleg Gauss-Newton'}
                    tic
                    [xsol,flag,logg] = doglet_gauss_newton(obj);
                    t = toc;
                    logg = logg.iter_dataprocessing;        
                case {'Steepest Descent'}
                    [xsol,flag,logg] = steepest_descent(obj);
            
                case {'internal'}
                    opt = optimoptions('lsqnonlin','SpecifyObjectiveGradient',true,'Display','off');
                    [xsol,~,~,flag] = lsqnonlin(@(x)myfun(x,obj.r,obj.drdx), obj.x0, [],[],opt);
                     logg            = nlsOption;
            end
        end      
        % check the dimension of lower and upper bounds
        function bool = check_bound_dim(obj)
            Nx = numel(obj.x0);
            if isempty(obj.bound)
                warning('no lower and upper bounds')
            else
                sz = size(obj.bound);
                if sz(1) == Nx && sz(2) == 2
                    bool = true;
                else
                    bool = false;
                    warning('bound dimension is wrong')
                end
            end
        end
    end
end

function [F,J] = myfun(x,r,dr)
F = r(x);     % Objective function values at x
if nargout > 1   % Two output arguments
   J = dr(x);   % Jacobian of the function evaluated at x
end
end
