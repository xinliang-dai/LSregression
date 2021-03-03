function [d,y,f,dfdx,x0] = dim_1_case
% initial one-dimension test
    d = linspace(0,3)';
    y = exp(-1.3*d) + 0.05*randn(size(d));
    x0= 5;
    f    = @(r)exp(-d*r)-y;
    dfdx = @(r)-d .* exp(-d*r);
end