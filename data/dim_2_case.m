function [d,y,f,dfdx,x0,xsol] = dim_2_case
% initial one-dimension test
    d    = linspace(0,3)';
    xsol = [-1.3; 0.9];
%     xsol = randn([2,1]);
    sigma = norm(xsol,inf);
    y    = exp(-1.3*d) - sin(0.9*d) ;%+ 0.1*randn(size(d));
    
    x0   = vary_initial_point(xsol,sigma);

%     x0   = [-1.3; 0.9];
    f    = @(x)exp(d*x(1)) - sin(d*x(2))-y;
    dfdx = @(x)[d .* exp(d*x(1)), -d.* cos(d*x(2))];
end