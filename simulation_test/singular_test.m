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
[data,r,drdx,cost] = singular_Powell;


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


data.post_dataprocessing(test_model.r);
% iter Information processing
%% results
data.post_dataprocessing(test_model.r);


%%
nls_method  = 'preconditioned Levenberg-Marquardt';
cg_method  = 'CG-Steihaug';
test_model.options = nlsOption(iter_max, tol, nls_method, cg_method);
[data.xsol, data.logg, flag] = test_model.solve_nls;                    

data.post_dataprocessing(test_model.r);

%%
nls_method  = 'Gauss-Newton';
cg_method  = 'preconditioned CG-Steihaug';
test_model.options = nlsOption(iter_max, tol, nls_method, cg_method);
[data.xsol, data.logg, flag] = test_model.solve_nls;                    
track = [data.x0,data.logg.xk];
plot(track(1,:),track(2,:))

data.post_dataprocessing(test_model.r);




