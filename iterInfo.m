classdef iterInfo
    % record information in iterations
    properties
        iter           int8 {mustBeNumeric}
        xk             double {mustBeNumeric}
        fval           double {mustBeNumeric}
        dfval          double {mustBeNumeric}
        rel_steplength double {mustBeNumeric}
    end
    
    methods
        % constructor with maximal iteration and dimension of variables
        function obj = iterInfo(iter_max,Nx)
            if nargin>1
                obj.xk = zeros(Nx,iter_max);
            end
            % initial recording
            obj.iter           = zeros(iter_max,1);
            obj.fval           = zeros(iter_max,1);
            obj.dfval          = zeros(iter_max,1);
            obj.rel_steplength = zeros(iter_max,1);
        end
        
        % post-dataprocessing: reduce zeros column
        function obj = post_dataprocessing(obj)
            if isempty(obj.iter)
                warning('iterative record error')
            else
                % reduce dimension                
                idx         = find(obj.iter);
                obj.iter    = obj.iter(idx);
                obj.xk      = obj.xk(:,idx);
                obj.fval    = obj.fval(idx);
                obj.dfval   = obj.dfval(idx);
                obj.rel_steplength = obj.rel_steplength(idx);
            end
        end

    end
end

