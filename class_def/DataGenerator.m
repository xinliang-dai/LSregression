classdef DataGenerator
    %DATAGENERATOR Summary of this class goes here
    %   Detailed explanation goes here    
    properties
        d               {mustBeNumeric}
        y               {mustBeNumeric}
        x0              {mustBeNumeric}
        xref            {mustBeNumeric}
        xsol            {mustBeNumeric}
        logg   iterInfo        
    end
    
    methods
        function post_dataprocessing(obj,r)
            %% data processing of simlation results
   
            % compare result with ref solution
            e = norm(obj.xref-obj.xsol,2)
            % plot iter information
            plot_results(obj,r)
        end
    end
end


function plot_results(obj,r)

xref  = obj.xref;
logg  = obj.logg;
d     = obj.d;
y     = obj.y;
Niter = numel(logg.iter);
dx    = zeros(Niter,1);
for i = 1:Niter
    dx(i) = norm(logg.xk(:,i)-xref,2);
    if dx(i) ==0
        dx(i) = eps;
    end
end

x_bound = [1, logg.iter(end)+1];

figure('Name','Iter Plotting')
subplot(3,2,1)
semilogy(logg.iter, dx, '^-')
xlabel('$\mathrm{Iteration}$','fontsize',12,'interpreter','latex')
ylabel('$||x^*-x_k||_2$','fontsize',12,'interpreter','latex')
xlim(x_bound)
ylim(plot_y_limit(dx))
grid on 

subplot(3,2,2)
semilogy(logg.iter, logg.fval, '^-')
xlabel('$\mathrm{Iteration}$','fontsize',12,'interpreter','latex')
ylabel('$f(x_k)$','fontsize',12,'interpreter','latex')
xlim(x_bound)
ylim(plot_y_limit(logg.fval))
grid on 

subplot(3,2,3)
logg.grad_norm(logg.grad_norm==0) = eps;
semilogy(logg.iter, logg.grad_norm, '^-')
xlabel('$\mathrm{Iteration}$','fontsize',12,'interpreter','latex')
ylabel('$||\nabla f(x_k)||_2$','fontsize',12,'interpreter','latex')
xlim(x_bound)
ylim(plot_y_limit(logg.grad_norm))

grid on 

subplot(3,2,4)
logg.dfval(logg.dfval==0) = eps;
semilogy(logg.iter, logg.dfval, '^-')
xlabel('$\mathrm{Iteration}$','fontsize',12,'interpreter','latex')
ylabel('$f(x_k)-f(x_k+p_k)$','fontsize',12,'interpreter','latex')
xlim(x_bound)
ylim(plot_y_limit(logg.dfval))

grid on 


subplot(3,2,5)
semilogy(logg.iter, logg.delta, '^-')
xlabel('$\mathrm{Iteration}$','fontsize',12,'interpreter','latex')
ylabel('$\Delta$','fontsize',12,'interpreter','latex')
xlim(x_bound)
ylim(plot_y_limit(logg.delta))
grid on 

subplot(3,2,6)
logg.pk_norm(logg.pk_norm==0) = eps;
semilogy(logg.iter, logg.pk_norm, '^-')
xlabel('$\mathrm{Iteration}$','fontsize',12,'interpreter','latex')
ylabel('$||p_k||_2$','fontsize',12,'interpreter','latex')
xlim(x_bound)
ylim(plot_y_limit(logg.pk_norm))
grid on 

y_estimate = r(logg.xk(:,end))+y;
figure('Name','Best fit for measured data')
plot(d,y,'ko',d,y_estimate,'b-')
legend('Data','Best fit')
xlabel('t')
ylabel('data')
% ylim(plot_y_limit(y_estimate))
grid on 

end

function bound =  plot_y_limit(y)
    y_min = min(y)/10;
    y_max = max(y)*10;
    if y_min == 0
        y_min = eps;
        if y_max == 0
            y_max = 1;
        end
    end
    bound = [y_min, y_max];
end
