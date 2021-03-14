function [data,r,drdx] = dim_4_case
    % initial object of data
    data      = DataGenerator;
    % d vector
    t         = linspace(0,3)';
    data.t    = t;
    data.xref = [-1;-2;1;-1];
    sigma     = norm(data.xref,inf);
    % measured data with noise
    data.y    = data.xref(3)*exp(data.xref(1).*t)+data.xref(4)*exp(data.xref(2).*t);%+0.5*randn(size(d));
    % varying initial point
    data.x0   = [-4 -5 4 -4]';
%     data.x0   = vary_initial_point(data.xref,sigma);
    % residual vector-function
    r         = @(x) x(3)*exp(x(1)*t)+x(4)*exp(x(2)*t) -data.y;
    % first order derivative of residual - Jacobian matrix
    drdx      = @(x) [x(3)*t.*exp(x(1)*t), x(4)*t.*exp(x(2)*t),  exp(x(1)*t),exp(x(2)*t)];
end