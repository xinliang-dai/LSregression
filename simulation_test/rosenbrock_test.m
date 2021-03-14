%% add path
clc
clear
close all
addpath(genpath('../data/'));
addpath(genpath('../class_def/'));
addpath(genpath('../nls_solvers/'));
% create data 
% [data,r,drdx] = dim_1_case;
% [data,r,drdx] = dim_2_case;
% [data,r,drdx] = dim_2_quadratic_case;
[data,r,drdx,cost] = rosenbrock;

test_model      = nlsProblem(r, drdx, data.x0);

% setting options for solving nls problem
tol        = 10^(-10);
iter_max   = 1000;
nls_method = 'Gauss-Newton';
% nls_method = 'Levenberg-Marquardt';
% cg_method  = 'without CG';
% cg_method  = 'standard CG';
% cg_method  = 'preconditioned CG';
% cg_method  = 'CG-Steihaug';
% cg_method  = 'preconditioned CG-Steihaug';
% nls_method  = 'Levenberg-Marquardt';
% nls_method  = 'preconditioned Levenberg-Marquardt';
% cg_method  = 'without CG';

x = linspace(-3,3); y = linspace(-3,5);
[xx,yy] = meshgrid(x,y); 
ff = log(cost(xx,yy));
levels = 1:8;
LW = 'linewidth'; FS = 'fontsize'; MS = 'markersize';
figure, contour(x,y,ff,levels,LW,1.2), colormap(jet) 
colorbar
hold on
% axis([-1.5 1.5 -1 3]), axis square, hold on

% initial test model
test_model      = nlsProblem(r, drdx, data.x0);

% setting options for solving nls problem
tol        = 10^(-10);
iter_max   = 1000;
% nls_method = 'Gauss-Newton';
% nls_method = 'Levenberg-Marquardt';
% cg_method  = 'without CG';
% cg_method  = 'standard CG';
% cg_method  = 'preconditioned CG';
% cg_method  = 'CG-Steihaug';
% cg_method  = 'preconditioned CG-Steihaug';
% nls_method  = 'Levenberg-Marquardt';
% nls_method  = 'preconditioned Levenberg-Marquardt';
% cg_method  = 'without CG';
% cg_method  = 'CG-Steihaug';
nls_method  = 'Dogleg Gauss-Newton';
cg_method  = 'preconditioned CG';

%% run internal nls solver
test_model.options = nlsOption(iter_max, tol, nls_method, cg_method);
[data.xsol, data.logg, flag] = test_model.solve_nls;            
track = [data.x0,data.logg.xk];
plot(track(1,:),track(2,:))

data.post_dataprocessing(test_model.r);
% iter Information processing
%% results
data.post_dataprocessing(test_model.r);


%%
nls_method  = 'preconditioned Levenberg-Marquardt';
cg_method  = 'CG-Steihaug';
test_model.options = nlsOption(iter_max, tol, nls_method, cg_method);
[data.xsol, data.logg, flag] = test_model.solve_nls;                    
track = [data.x0,data.logg.xk];
plot(track(1,:),track(2,:))
data.post_dataprocessing(test_model.r);

%%
nls_method  = 'Gauss-Newton';
cg_method  = 'standard CG';
test_model.options = nlsOption(iter_max, tol, nls_method, cg_method);
[data.xsol, data.logg, flag] = test_model.solve_nls;                    
track = [data.x0,data.logg.xk];
plot(track(1,:),track(2,:))

data.post_dataprocessing(test_model.r);




