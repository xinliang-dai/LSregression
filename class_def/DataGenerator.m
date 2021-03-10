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
            % iter Information processing
            if ~isnan(obj.logg.iter) 
                obj.logg = obj.logg.iter_dataprocessing;
            end
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
end

x_bound = [1, logg.iter(end)+1];

figure('Name','Iter Plotting')
subplot(3,2,1)
semilogy(logg.iter, dx)
xlabel('$\mathrm{Iteration}$','fontsize',12,'interpreter','latex')
ylabel('$||x^*-x^k||_2$','fontsize',12,'interpreter','latex')
xlim(x_bound)
ylim(plot_y_limit(dx))
grid on 

subplot(3,2,2)
semilogy(logg.iter, logg.fval)
xlabel('$\mathrm{Iteration}$','fontsize',12,'interpreter','latex')
ylabel('$f(x^k)$','fontsize',12,'interpreter','latex')
xlim(x_bound)
ylim(plot_y_limit(logg.fval))
grid on 

subplot(3,2,3)
semilogy(logg.iter, logg.grad_norm)
xlabel('$\mathrm{Iteration}$','fontsize',12,'interpreter','latex')
ylabel('$||\nabla f(x^k)||_2$','fontsize',12,'interpreter','latex')
ylim(plot_y_limit(logg.grad_norm))

grid on 

subplot(3,2,4)
semilogy(logg.iter, logg.dfval)
xlabel('$\mathrm{Iteration}$','fontsize',12,'interpreter','latex')
ylabel('$m_k(x_k)-m_k(x_k+p_k)$','fontsize',12,'interpreter','latex')
ylim(plot_y_limit(logg.dfval))

grid on 


subplot(3,2,5)
semilogy(logg.iter, logg.delta)
xlabel('$\mathrm{Iteration}$','fontsize',12,'interpreter','latex')
ylabel('$\Delta$','fontsize',12,'interpreter','latex')
xlim(x_bound)
ylim(plot_y_limit(logg.delta))
grid on 

subplot(3,2,6)
semilogy(logg.iter, logg.pk_norm)
xlabel('$\mathrm{Iteration}$','fontsize',12,'interpreter','latex')
ylabel('$||p||_2$','fontsize',12,'interpreter','latex')
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
    if y_min == 0
        bound = [eps,  max(y)*10];
    else
        bound = [min(y)/10,  max(y)*10];
    end
end