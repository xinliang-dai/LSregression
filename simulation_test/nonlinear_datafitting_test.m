%% add path
clc
clear
close all
addpath(genpath('../data/'));
addpath(genpath('../class_def/'));
addpath(genpath('../nls_solvers/'));
%% create data 
% [data,r,drdx] = dim_1_case;
% [data,r,drdx] = dim_2_case;
% [data,r,drdx] = dim_2_quadratic_case;
[data,r,drdx] = dim_4_case(1,0);


%% initial math model for nonlinear least-squares problem
% initial test model
test_model      = nlsProblem(r, drdx, data.x0);
data.x0 = [-1.09135605225178
1.7057910722914
-3.63137589048434
-0.411043248578462]';
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
% iter Information processing
%% results
data.post_dataprocessing(test_model.r);
% 
% %%
% nls_method  = 'preconditioned Levenberg-Marquardt';
% cg_method  = 'CG-Steihaug';
% test_model.options = nlsOption(iter_max, tol, nls_method, cg_method);
% [data.xsol, data.logg, flag] = test_model.solve_nls;                    
% data.post_dataprocessing(test_model.r);
% 
% %%
% nls_method  = 'Gauss-Newton';
% cg_method  = 'preconditioned CG-Steihaug';
% test_model.options = nlsOption(iter_max, tol, nls_method, cg_method);
% [data.xsol, data.logg, flag] = test_model.solve_nls;                    
% data.post_dataprocessing(test_model.r);
%%
nls_method  = 'internal';
cg_method  = 'preconditioned CG-Steihaug';
test_model.options = nlsOption(iter_max, tol, nls_method, cg_method);
[data.xsol, ~, flag] = test_model.solve_nls;                    
% data.post_dataprocessing(test_model.r);
e      = norm(data.xref-data.xsol,2)
