%% create data 
% [d,y,f,dfdx,x0] = dim_1_case;
[d,y,f,dfdx,x0] = dim_2_case;
%% initial obj
test      = nlsProblem;
test.x0   = x0;
test.f    = f;
test.dfdx = dfdx;
% setting options for solving nls problem
tol       = 10^(-10);
iter_max  = 40;
test.options = nlsOption(iter_max,tol,'class gauss-newton');
%% run internal nls solver
[xsol_gn, logg, flag] = test.solve_nls;
xsol_gn
%% plot result
plot(d,y,'ko',d,test.f(xsol_gn)+y,'b-')
legend('Data','Best fit')
xlabel('t')
ylabel('exp(-tx)')