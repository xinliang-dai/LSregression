classdef nlsOption
    %NLSOPTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        iter_max  int8   {mustBePositive} 
        tol       double {mustBePositive} 
        method         char 
        solver         char  
    end
    
    methods
        % constructor
        function obj = nlsOption(iter_max,tol,method,solver)
            
            if nargin > 0
                obj.iter_max = iter_max;
                if nargin > 1
                    obj.tol      = tol;
                    if nargin >2
                        obj.method   = method;
                        if nargin > 3
                             obj.solver   = solver;
                        end
                    end
                end
            end

        end
        
%         function outputArg = method1(obj,inputArg)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
%         end
    end
end

