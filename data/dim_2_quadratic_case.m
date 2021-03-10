function [d,y,r,drdx,x0,xval] = dim_2_quadratic_case
% initial one-dimension test
    xval = [1 ;3];
    radiu = 1;
    sigma = 5;
    d = linspace(0,3)';
    y = xval(1)*d+xval(2)*d.^2 + 0.05*randn(size(d));
    x0= vary_initial_point(xval,sigma);
    r    = @(x) x(1)*d+x(2)*d.^2 -y;
    drdx = @(x) [d,d.^2];
end