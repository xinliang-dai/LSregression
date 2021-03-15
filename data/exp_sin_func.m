function [data,r,drdx,cost] = exp_sin_func
    % initial object of data
    data      = DataGenerator;
    % d vector
%     t         = linspace(0,5)';
%     data.t    = t;
    data.xref = [-1.3; 1.3];
    sigma     = norm(data.xref,inf);
    % measured data with noise
%     data.y    = exp(-1.3*t) .* sin(1.3*t) +0.01 * sigma * randn(size(t));
    % varying initial point
    data.x0   = [1;1];% vary_initial_point(data.xref,sigma);
    % residual vector-function
    r         = @(x)exp(x(1)) + sin(x(2));
    cost      = @(x,y)exp(x).^2 .* y.^2/2;
    % first order derivative of residual - Jacobian matrix
    drdx      = @(x)[exp(x(1)) .* x(2), exp(x(1))];
end