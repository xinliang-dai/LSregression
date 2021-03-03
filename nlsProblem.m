classdef nlsProblem
    %   define standard nonlinear problem
    %   the nonlinear least squares objective function F(x) := 0.5*(f(x)'*f(x)).
    properties
        % residual
        f 
        % jacobian
        dfdx 
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
                case {'internal solver'}
                    [xsol,~,~,flag] = lsqnonlin(obj.f, obj.x0);
                    logg            = nan;
                case {'classic gauss-newton'}
                    [xsol,flag,logg] = basic_gauss_newton(obj);
                    logg = logg.post_dataprocessing;
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

