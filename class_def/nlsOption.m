classdef nlsOption
    %NLSOPTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        iter_max           {mustBePositive} 
        tol         double {mustBePositive} 
        cg_method   char 
        nls_method  char  
    end
    
    methods
        % constructor
        function obj = nlsOption(iter_max,tol,nls_method,cg_method)            
            if nargin > 0 
                obj.iter_max   = iter_max; 
            else
                obj.iter_max   = 200;
            end                
            if nargin > 1 
                obj.tol        = tol;
            else
                obj.tol        = 10^-6;
            end
            if nargin >2
                obj.nls_method = nls_method;
            else
                obj.nls_method = 'Gauss-Newton';                        
            end
            if nargin >3
                obj.cg_method  = cg_method;
            else
                obj.cg_method  = 'preconditioned CG-Steihaug';
            end
        end
    end
end

