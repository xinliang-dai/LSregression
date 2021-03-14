function [data,r,drdx] = dim_1_case
    % initial object of data
    data      = DataGenerator;
    % d vector
    t         = linspace(0,5)';
    data.d    = t;
    data.xref = -1.3;
    sigma     = abs(data.xref)/2;
    % measured data with noise
    data.y    = exp(data.xref*t) + 0.05*randn(size(t));
    % varying initial point
    data.x0   = vary_initial_point(data.xref,sigma);
    % residual vector-function
    r         = @(x)exp(t*x)-data.y;
    % first order derivative of residual - Jacobian matrix
    drdx      = @(x)t .* exp(t*x);
end