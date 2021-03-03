%% create data 
d = linspace(0,3)';
y = exp(-1.3*d) + 0.05*randn(size(d));
% y = 1.3*d + 0.05*randn(size(d));
%% initial obj
test      = nlsProblem;
test.x0   = 5;
test.f    =  @(r)exp(-d*r)-y;
test.dfdx =  @(r)-exp(-d*r);

tol       = 10^(-3);
iter_max  = 40;
test.options = nlsOption(iter_max,tol,'classic gauss-newton');
% test.f = @(x) (d*x-y);
% test.dfdx =  @(x) d;
%% run internal nls solver
% xsol = test.internal_solver
% [xsol_gn,flag, logg] = basic_gauss_newton(test);
[xsol_gn, logg, flag] = test.solve_nls;
% post date-processing
%% plot result
plot(d,y,'ko',d,d*xsol_gn,'b-')
legend('Data','Best fit')
xlabel('t')
ylabel('exp(-tx)')