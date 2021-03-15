function [data,r,drdx] = dim_4_case(mean,theta)
    % initial object of data
    data      = DataGenerator;
    % d vector
    t         = linspace(0,3)';
    data.t    = t;
    data.xref = [-1;-2;1;-1];
    sigma     = norm(data.xref,inf);
    % measured data with noise
    data.y    = data.xref(3)*exp(data.xref(1).*t)+data.xref(4)*exp(data.xref(2).*t)+theta*randn(size(t));
    % varying initial point
    Nx        = numel(data.xref);
    data.x0   = data.xref + vary_initial(Nx,mean,sigma);
%     data.x0   = vary_initial_point(data.xref,sigma);
    % residual vector-function
    r         = @(x) x(3)*exp(x(1)*t)+x(4)* exp(x(2)*t) -data.y;
    % first order derivative of residual - Jacobian matrix
    drdx      = @(x) [x(3)*t.*exp(x(1)*t), x(4)*t.*exp(x(2)*t),  exp(x(1)*t),exp(x(2)*t)];
end

function x0 = vary_initial(Nx,radiu,sigma)
        random_vector = (2*rand(Nx,1)-1)*radiu;
        x0=random('Normal',random_vector,sigma*ones(Nx,1)); 
end