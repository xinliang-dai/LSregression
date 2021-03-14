function [data,r,drdx] = dim_2_quadratic_case
    % initial object of data
    data      = DataGenerator;
    % d vector
    t         = linspace(0,5)';
    data.d    = t;
    data.xref = [1 ;3];
    sigma     = norm(data.xref,inf);
    % measured data with noise
    data.y    = data.xref(1)*t+data.xref(2)*t.^2 + 1*randn(size(t));
    % varying initial point
    data.x0   = vary_initial_point(data.xref,sigma);
    % residual vector-function
    r         = @(x) x(1)*t+x(2)*t.^2 -data.y;
    % first order derivative of residual - Jacobian matrix
    drdx      = @(x) [t,t.^2];
end