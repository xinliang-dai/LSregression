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
        function [xsol, logg, flag]= solve_nls(obj)
            % switch from different methods
            if isempty(obj.options)
                obj.options = nlsOption;
            end
            switch obj.options.nls_method
                case {'Gauss-Newton'}
                    [xsol,flag,logg] = pcg_steihaug_gauss_newton(obj);
                    logg = logg.iter_dataprocessing;
                case {'Levenberg-Marquardt'}
                    [xsol,flag,logg] = standart_levenberg_marquardt(obj);
                    logg = logg.iter_dataprocessing;                    
            end
%                 case {'internal'}
%                     [xsol,~,~,flag] = lsqnonlin(obj.r, obj.x0);
%                     logg            = nan;
%                 case {'class gauss-newton'}
%                     [xsol,flag,logg] = basic_gauss_newton(obj);
% %                     logg = logg.iter_dataprocessing;
%                 case ('CG gauss-newton')
%                     [xsol,flag,logg] = cg_gauss_newton(obj);
% %                     logg = logg.iter_dataprocessing;  
%                 case ('PCG gauss-newton')
%                     [xsol,flag,logg] = pcg_gauss_newton(obj);
% %                     logg = logg.iter_dataprocessing;   
%                 case ('CG-Steihaug gauss-newton')
%                     [xsol,flag,logg] = cg_steihaug_gauss_newton(obj);
% %                     logg = logg.iter_dataprocessing;
%                 case ('PCG-Steihaug gauss-newton')
%                     [xsol,flag,logg] = pcg_steihaug_gauss_newton(obj);
% %                     logg = logg.iter_dataprocessing;
%                             
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

% function [f,df]

