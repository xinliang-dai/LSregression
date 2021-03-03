classdef nlsOption
    %NLSOPTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        iter_max (1,1) int8   {mustBePositive} = 50;
        tol      (1,1) double {mustBePositive} = 10^(-3);
        method         char = 'classic gauss-newton'
        solver         char = 'internal'  
    end
    
    methods
        % constructor
        function obj = nlsOption(iter_max,tol,method,solver)
            if nargin>4
                obj.solver   = solver;
            elseif nargin>3
                obj.method   = method;
            elseif nargin>2
                obj.tol      = tol;
            elseif nargin>1
                obj.iter_max = iter_max;
            end          
        end
        
%         function outputArg = method1(obj,inputArg)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
%         end
    end
end

