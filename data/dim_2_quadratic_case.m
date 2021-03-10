function [data,r,drdx] = dim_2_quadratic_case
    % initial object of data
    data      = DataGenerator;
    % d vector
    d         = linspace(0,5)';
    data.d    = d;
    data.xref = [1 ;3];
    sigma     = norm(data.xref,inf);
    % measured data with noise
    data.y    = data.xref(1)*d+data.xref(2)*d.^2 + 1*randn(size(d));
    % varying initial point
    data.x0   = vary_initial_point(data.xref,sigma);
    % residual vector-function
    r         = @(x) x(1)*d+x(2)*d.^2 -data.y;
    % first order derivative of residual - Jacobian matrix
    drdx      = @(x) [d,d.^2];
end