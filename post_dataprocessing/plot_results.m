function plot_results(xref,logg,f,d,y)
Niter = numel(logg.iter);
dx    = zeros(Niter,1);
for i = 1:Niter
    dx(i) = norm(logg.xk(:,i)-xref,2);
end

figure('Name','Iter Plotting')
subplot(3,2,1)
semilogy(logg.iter, dx)
xlabel('$\mathrm{Iteration}$','fontsize',12,'interpreter','latex')
ylabel('$||x^*-x^k||_2$','fontsize',12,'interpreter','latex')
grid on 

subplot(3,2,2)
semilogy(logg.iter, logg.fval)
xlabel('$\mathrm{Iteration}$','fontsize',12,'interpreter','latex')
ylabel('$f(x^k)$','fontsize',12,'interpreter','latex')

grid on 

subplot(3,2,3)
semilogy(logg.iter, logg.dfval)
xlabel('$\mathrm{Iteration}$','fontsize',12,'interpreter','latex')
ylabel('$\nabla f(x^k)$','fontsize',12,'interpreter','latex')

grid on 

subplot(3,2,4)
semilogy(logg.iter, logg.pk_norm)
xlabel('$\mathrm{Iteration}$','fontsize',12,'interpreter','latex')
ylabel('$\frac{||p||_2}{||x^k||_2}$','fontsize',12,'interpreter','latex')

grid on 


subplot(3,2,5)
semilogy(logg.iter, logg.delta)
xlabel('$\mathrm{Iteration}$','fontsize',12,'interpreter','latex')
ylabel('$\frac{||p||_2}{||x^k||_2}$','fontsize',12,'interpreter','latex')

grid on 


figure('Name','Best fit for measured data')
plot(d,y,'ko',d,f(logg.xk(:,end))+y,'b-')
legend('Data','Best fit')
xlabel('t')
ylabel('data')
end