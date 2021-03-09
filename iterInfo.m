classdef iterInfo
    % record information in iterations
    properties
        iter           int8 {mustBeNumeric}
        xk             double {mustBeNumeric}
        fval           double {mustBeNumeric}
        dfval          double {mustBeNumeric}
        pk             double {mustBeNumeric}
        delta          double {mustBeNumeric}
        xk_norm        double {mustBeNumeric}
        pk_norm        double {mustBeNumeric}        
    end
    
    methods
        % constructor with maximal iteration and dimension of variables
        function obj = iterInfo(iter_max,Nx)
            if nargin>1
                obj.xk = zeros(Nx,iter_max);
                obj.pk = zeros(Nx,iter_max);
            end
            % initial recording
            obj.iter           = zeros(iter_max,1);
            obj.fval           = zeros(iter_max,1);
            obj.dfval          = zeros(iter_max,1);
            obj.delta          = zeros(iter_max,1);
        end
        
        % post-dataprocessing: reduce zeros column
        function obj = post_dataprocessing(obj)
            if isempty(obj.iter)
                warning('iterative record error')
            else
                % reduce dimension                
                idx         = find(obj.iter);
                obj.iter    = obj.iter(idx);
                obj.fval    = obj.fval(idx);
                obj.dfval   = obj.dfval(idx);
                if ~isempty(obj.delta)
                    obj.delta   = obj.delta(idx);
                end
                if ~isempty(obj.xk) && ~isempty(obj.pk)
                    obj.xk      = obj.xk(:,idx);
                    obj.pk      = obj.pk(:,idx);
                    % norm of xk and pk
                    for i = 1:numel(idx)
                        obj.xk_norm(i) = obj.xk(:,i)'* obj.xk(:,i);
                        obj.pk_norm(i) = obj.pk(:,i)'* obj.pk(:,i);
                    end
                end
            end
        end
    end
end

