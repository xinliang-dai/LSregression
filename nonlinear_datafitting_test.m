%% add path
clear;
close all;
clc
addpath(genpath('../nonlinear_LSsolvers/'));
%% create data 
% [d,y,r,drdx,x0,xval] = dim_1_case;
[d,y,r,drdx,x0,xval] = dim_2_case;
% [d,y,r,drdx,x0,xval] = dim_2_quadratic_case;
%% initial obj
test      = nlsProblem;
test.x0   = x0;
test.r    = r;
test.drdx = drdx;
% setting options for solving nls problem
tol       = 10^(-6);
iter_max  = 100;
test.options = nlsOption(iter_max,tol,'CG-Steihaug gauss-newton');
%% run internal nls solver
[xsol_gn, logg, flag] = test.solve_nls;
logg = logg.post_dataprocessing;                       

%% plot result
compare_results(xval,xsol_gn)
plot_results(xval,logg,r,d,y)
