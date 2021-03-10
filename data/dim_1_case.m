function [data,r,drdx] = dim_1_case
% % initial one-dimension test
%     xval = -1.3;
%     radiu = 1;
%     sigma = 1;
%     d = linspace(0,3)';
%     y = exp(xval*d) + 0.05*randn(size(d));
%     x0= vary_initial_point(xval,sigma);
%     f    = @(x)exp(d*x)-y;
%     dfdx = @(x)d .* exp(d*x);
%     
    % initial object of data
    data      = DataGenerator;
    % d vector
    d         = linspace(0,5)';
    data.d    = d;
    data.xref = -1.3;
    sigma     = abs(data.xref)/2;
    % measured data with noise
    data.y    = exp(data.xref*d) + 0.05*randn(size(d));
    % varying initial point
    data.x0   = vary_initial_point(data.xref,sigma);
    % residual vector-function
    r         = @(x)exp(d*x)-data.y;
    % first order derivative of residual - Jacobian matrix
    drdx      = @(x)d .* exp(d*x);
end