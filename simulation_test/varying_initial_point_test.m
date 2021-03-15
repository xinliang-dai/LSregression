%% add path
clc
clear
close all
addpath(genpath('../data/'));
addpath(genpath('../class_def/'));
addpath(genpath('../nls_solvers/'));
%% create data 
tol             = 10^(-6);
iter_max        = 200;
n               = 100;
iter_GND        = zeros(n,1);
t_GND           = zeros(n,1);
success_GND     = zeros(n,1);
e_GND           = zeros(n,1);
iter_LM        = zeros(n,1);
t_LM           = zeros(n,1);
success_LM     = zeros(n,1);
e_LM           = zeros(n,1);

iter_GN        = zeros(n,1);
t_GN           = zeros(n,1);
success_GN     = zeros(n,1);
e_GN           = zeros(n,1);
x0             = zeros(4,n);

dx_norminf     = zeros(n,1);
dx_norm2       = zeros(n,1);

theta = 0;

xref           = [-2;-1;-1;1];
e_internal     = zeros(n,1);
t_internal     = zeros(n,1);
for i = 5:5:100
    for j = 1:n
        [data,r,drdx]   = dim_2_case(i,theta);
        test_model      = nlsProblem(r, drdx, data.x0);
%         x0(:,j)         = data.x0;
        dx_norminf(j)   = min(norm(data.x0-data.xref,inf));
        dx_norm2(j)     = min(norm(data.x0-data.xref,2));
        nls_method      = 'internal';
        cg_method       = 'standard CG';
        test_model.options = nlsOption(iter_max, tol, nls_method, cg_method);
        [data.xsol, ~, flag] = test_model.solve_nls;                    
        e_internal(j)        = norm(data.xsol-data.xref,2);
        if e_internal(j)<=1e-5
            success_GND(j)=true;
        else
            success_GND(j)=false;
        end
        
        %% Gauss-Newton Dogleg
        nls_method      = 'Dogleg Gauss-Newton';
        cg_method       = 'standard CG';
        test_model.options = nlsOption(iter_max, tol, nls_method, cg_method);
        [data.xsol, data.logg, flag, t_GND(j)] = test_model.solve_nls;                    
        iter_GND(j)     = numel(data.logg.iter);
        e_GND(j)        = norm(data.xsol-data.xref,2);
        if e_GND(j)<=1e-5
            success_GND(j)=true;
        else
            success_GND(j)=false;
        end
        
        %% Levenberg-Marquardt
        nls_method      = 'Levenberg-Marquardt';
        cg_method       = 'CG-Steihaug';
        test_model.options = nlsOption(iter_max, tol, nls_method, cg_method);
        [data.xsol, data.logg, flag,  t_LM(j)] = test_model.solve_nls;                    
        iter_LM(j)      = numel(data.logg.iter);
        e_LM(j)         = norm(data.xsol-data.xref,2);
        if e_LM(j)<=1e-5
            success_LM(j)=true;
        else
            success_LM(j)=false;
        end
        %% Gauss-Newton preconditioned CG-Steihaug
        nls_method     = 'Gauss-Newton';
        cg_method      = 'CG-Steihaug';
        test_model.options = nlsOption(iter_max, tol, nls_method, cg_method);
        [data.xsol, data.logg, flag,  t_GN(j)] = test_model.solve_nls;
        iter_GN(j)     = numel(data.logg.iter);
        e_GN(j)        = norm(data.xsol-data.xref,2);
        if e_GN(j)<=1e-5
            success_GN(j)=true;
        else
            success_GN(j)=false;
        end
    end
    % save results to table
    file_name = "vary_initial_record\Results_GND_"+i+".csv";
    T_results = table(dx_norminf,dx_norm2,success_GND,e_GND, iter_GND,t_GND);
    writetable(T_results,file_name)
    % save results to table
    file_name = "vary_initial_record\Results_LM_"+i+".csv";
    T_results = table(dx_norminf,dx_norm2,success_LM,e_LM, iter_LM, t_LM);
    writetable(T_results,file_name)
    % save results to table
    file_name = "vary_initial_record\Results_GN_"+i+".csv";
    T_results = table(dx_norminf,dx_norm2,success_GN,e_GN, iter_GN, t_GN);
    writetable(T_results,file_name)    
    % save initial dual to table
    file_name = "vary_initial_record\x0_dis_"+i+".csv";
    T_lam0 = table(x0);
    writetable(T_lam0,file_name)
    % save initial dual to table
    file_name = "vary_initial_record\flag_"+i+".csv";
     T    = table(success_GND,success_GND,success_LM,success_GN);
    writetable(T,file_name)    
end