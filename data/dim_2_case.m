function [data,r,drdx] = dim_2_case
    % initial object of data
    data      = DataGenerator;
    % d vector
    t         = linspace(0,5)';
    data.d    = t;
    data.xref = [-1.3; 0.9];
    sigma     = norm(data.xref,inf);
    % measured data with noise
    data.y    = exp(-1.3*t) .* sin(0.9*t) +0.01 * sigma * randn(size(t));
    % varying initial point
    data.x0   = vary_initial_point(data.xref,sigma);
    % residual vector-function
    r         = @(x)exp(t*x(1)) .* sin(t*x(2))-data.y;
    % first order derivative of residual - Jacobian matrix
    drdx      = @(x)[t .* exp(t*x(1)) .* sin(t*x(2)), t.* cos(t*x(2)) .*exp(t*x(1))];
end