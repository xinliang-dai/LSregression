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
        bound(:,2) double    {mustBeNumeric}
        options    nlsOption = nlsOption;
    end
    
    methods
        % solve the nonlinear least-squares problem
        function [xsol, logg, flag]= solve_nls(obj)
            obj.check_bound_dim;
            % switch from different methods
            switch obj.options.method
                case {'internal'}
                    [xsol,~,~,flag] = lsqnonlin(obj.r, obj.x0);
                    logg            = nan;
                case {'class gauss-newton'}
                    [xsol,flag,logg] = basic_gauss_newton(obj);
                    logg = logg.post_dataprocessing;
                case ('CG gauss-newton')
                    [xsol,flag,logg] = cg_gauss_newton(obj);
                    logg = logg.post_dataprocessing;  
                case ('PCG gauss-newton')
                    [xsol,flag,logg] = pcg_gauss_newton(obj);
                    logg = logg.post_dataprocessing;   
                case ('CG-Steihaug gauss-newton')
                    [xsol,flag,logg] = cg_steihaug_gauss_newton(obj);
            end                
        end      
        
        % check the dimension of lower and upper bounds
        function bool = check_bound_dim(obj)
            Nx = numel(obj.x0);
            if isempty(obj.bound)
                warning('no lower and upper bounds')
            else
                if size(obj.bound,1) == Nx
                    bool = true;
                else
                    bool = false;
                    warning('bound dimension is wrong')
                end
            end
        end
    end
end

% function [f,df]

