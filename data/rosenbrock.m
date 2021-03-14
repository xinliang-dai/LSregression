function [data,r,drdx,cost] = rosenbrock
    % initial object of data
    data      = DataGenerator;
    % d vector
%     d         = linspace(0,5)';
%     data.d    = d;
    data.xref = [1; 1];
    sigma     = norm(data.xref,inf);
    % measured data with noise
%     data.y    = 10;
    % varying initial point
    data.x0   = [-1;4];%vary_initial_point(data.xref,sigma);
    % residual vector-function
    r         = @(x) [10*(x(2)-x(1).^2); 1-x(1)];
    cost      =  @(x,y) (1-x).^2/2 + 100*(y-x.^2).^2/2;
    % first order derivative of residual - Jacobian matrix
    drdx      = @(x)[-20*x(1),10;
                     -1, 0];
end