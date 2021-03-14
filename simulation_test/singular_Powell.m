function [data,r,drdx,cost] = singular_Powell
    % initial object of data
    data      = DataGenerator;
    % d vector
%     d         = linspace(0,5)';
%     data.d    = d;
    data.xref = [0; 0];
    sigma     = norm(data.xref,inf);
    % measured data with noise
%     data.y    = 10;
    % varying initial point
    data.x0   = [3;1];%vary_initial_point(data.xref,sigma);
    % residual vector-function
    r         = @(x) [x(1);10*x(1)/(x(1)+0.1)+2*x(2)^2];
    cost      = @(x) x(1)^2/2+ 500*x(1)^2/((x(1)+0.1)+2*x(2)^2)^2;
    % first order derivative of residual - Jacobian matrix
    drdx      = @(x)[1,0;
                     (x(1)+0.1)^-2, 4*x(2)];
end