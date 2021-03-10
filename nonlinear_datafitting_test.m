%% add path
clc
clear
close all
addpath class_def;
addpath nls_solvers;
addpath data;
%% create data 
% [d,y,r,drdx,x0,xval] = dim_1_case;
[data,r,drdx] = dim_2_case;
% [d,y,r,drdx,x0,xval] = dim_2_quadratic_case;

%% initial math model for nonlinear least-squares problem
% initial test model
test_model      = nlsProblem(r, drdx, data.x0);

% setting options for solving nls problem
tol       = 10^(-10);
iter_max  = 100;

%% run internal nls solver
% solver_name        = 'gauss-newton with conjugate gradient';
% solver_name        = 'gauss-newton with preconditioned conjugate gradient';
% solver_name        = 'gauss-newton with CG-Steihaug ';
solver_name        = 'gauss-newton with PCG-Steihaug ';
test_model.options = nlsOption(iter_max,tol,'PCG-Steihaug gauss-newton');
[data.xsol, data.logg, flag] = test_model.solve_nls;                    
           
%% results
data.post_dataprocessing(test_model.r);