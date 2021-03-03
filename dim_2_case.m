function [d,y,f,dfdx,x0] = dim_2_case
% initial one-dimension test
    d = linspace(0,3)';
    y = exp(-1.3*d) - sin(0.9*d) + 0.05*randn(size(d));
    x0= [1;0.7];
    f    = @(r)exp(-d*r(1))-sin(d*r(2))-y;
    dfdx = @(r)-d .* exp(-d*r(1))- d.* cos(d*r(2));
end