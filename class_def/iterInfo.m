classdef iterInfo
    % record information in iterations
    properties
        iter           int8   {mustBeNumeric}
        xk             double {mustBeNumeric}
        fval           double {mustBeNumeric}
        grad           double {mustBeNumeric}
        dfval          double {mustBeNumeric}
        dmval          double {mustBeNumeric} 
        pk             double {mustBeNumeric}
        xk_norm        double {mustBeNumeric}
        pk_norm        double {mustBeNumeric}     
        grad_norm      double {mustBeNumeric}     
        rho            double {mustBeNumeric}     
        mu             double {mustBeNumeric}     
        delta          double {mustBeNumeric} 
    end
    
    methods
        % constructor with maximal iteration and dimension of variables
        function obj = iterInfo(iter_max,Nx)
            if nargin>1
                obj.xk   = zeros(Nx,iter_max);
                obj.pk   = zeros(Nx,iter_max);
                obj.grad = zeros(Nx,iter_max);
            end
            % initial recording
            obj.iter           = zeros(iter_max,1);
            obj.fval           = zeros(iter_max,1);
            obj.dfval          = zeros(iter_max,1);
            obj.dmval          = zeros(iter_max,1);
            obj.delta          = zeros(iter_max,1);
            obj.rho            = zeros(iter_max,1);
            obj.mu             = zeros(iter_max,1);
        end
        
        % post-dataprocessing: reduce zeros column
        function obj = iter_dataprocessing(obj)
            if isempty(obj.iter)
                warning('iterative record error')
            else
                % reduce dimension                
                idx         = find(obj.iter);
                obj.iter    = obj.iter(idx);
                obj.fval    = obj.fval(idx);
                obj.dfval   = obj.dfval(idx);
                obj.delta   = obj.delta(idx);
                obj.mu      = obj.mu(idx);
                obj.rho     = obj.rho(idx);
                if ~isempty(obj.xk) && ~isempty(obj.pk) && ~isempty(obj.grad)
                    obj.xk      = obj.xk(:,idx);
                    obj.pk      = obj.pk(:,idx);
                    obj.grad    = obj.grad(:,idx);
                    % norm of xk and pk
                    for i = 1:numel(idx)
                        obj.xk_norm(i)   = sqrt(obj.xk(:,i)'* obj.xk(:,i));
                        obj.pk_norm(i)   = sqrt(obj.pk(:,i)'* obj.pk(:,i));
                        obj.grad_norm(i) = sqrt(obj.grad(:,i)'* obj.grad(:,i));
                    end
                end
            end
        end
    end
end

