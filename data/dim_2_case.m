function [data,r,drdx] = dim_2_case
    % initial object of data
    data      = DataGenerator;
    % d vector
    d         = linspace(0,5)';
    data.d    = d;
    data.xref = [-1.3; 0.9];
    sigma     = norm(data.xref,inf);
    % measured data with noise
    data.y    = exp(-1.3*d) .* sin(0.9*d) +0.01 * sigma * randn(size(d));
    % varying initial point
    data.x0   = vary_initial_point(data.xref,sigma);
    % residual vector-function
    r         = @(x)exp(d*x(1)) .* sin(d*x(2))-data.y;
    % first order derivative of residual - Jacobian matrix
    drdx      = @(x)[d .* exp(d*x(1)) .* sin(d*x(2)), d.* cos(d*x(2)) .*exp(d*x(1))];
end