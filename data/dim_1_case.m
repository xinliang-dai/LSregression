function [d,y,f,dfdx,x0,xval] = dim_1_case
% initial one-dimension test
    xval = -1.3;
    radiu = 1;
    sigma = 1;
    d = linspace(0,3)';
    y = exp(xval*d) + 0.05*randn(size(d));
    x0= vary_initial_point(xval,sigma);
    f    = @(x)exp(d*x)-y;
    dfdx = @(x)d .* exp(d*x);
end