function [data,r,drdx] = dim_2_case(mean,theta)
    % initial object of data
    data      = DataGenerator;
    % d vector
    t         = linspace(0,5)';
    data.t    = t;
    data.xref = [-1.3; 1.3];
    sigma     = norm(data.xref,inf);
    % measured data with noise
    data.y    = exp(-1.3*t) .* sin(1.3*t) + theta * randn(size(t));
    % varying initial point
    Nx        = numel(data.xref);
    data.x0   = data.xref + vary_initial(Nx,mean,sigma);
    % residual vector-function
    r         = @(x)exp(t*x(1)) .* sin(t*x(2))-data.y;
    % first order derivative of residual - Jacobian matrix
    drdx      = @(x)[t .* exp(t*x(1)) .* sin(t*x(2)), t.* cos(t*x(2)) .*exp(t*x(1))];
end

function x0 = vary_initial(Nx,radiu,sigma)
        random_vector = (2*rand(Nx,1)-1)*radiu;
        x0=random('Normal',random_vector,sigma*ones(Nx,1)); 
end