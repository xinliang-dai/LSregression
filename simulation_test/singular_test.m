% %% add path
% clc
% clear
% close all
% addpath(genpath('../data/'));
% addpath(genpath('../class_def/'));
% addpath(genpath('../nls_solvers/'));
% % create data 
% % [data,r,drdx] = dim_1_case;
% % [data,r,drdx] = dim_2_case;
% % [data,r,drdx] = dim_2_quadratic_case;
% [data,r,drdx,cost] = singular_Powell;
% 
% 
% test_model      = nlsProblem(r, drdx, data.x0);
% 
% % setting options for solving nls problem
% tol        = 10^(-10);
% iter_max   = 1000;
% % nls_method = 'Gauss-Newton';
% % nls_method = 'Levenberg-Marquardt';
% % cg_method  = 'without CG';
% % cg_method  = 'standard CG';
% % cg_method  = 'preconditioned CG';
% % cg_method  = 'CG-Steihaug';
% % cg_method  = 'preconditioned CG-Steihaug';
% % nls_method  = 'Levenberg-Marquardt';
% % nls_method  = 'preconditioned Levenberg-Marquardt';
% % cg_method  = 'without CG';
% % cg_method  = 'CG-Steihaug';
% nls_method  = 'Dogleg Gauss-Newton';
% cg_method  = 'preconditioned CG';
% 
% %% run internal nls solver
% test_model.options = nlsOption(iter_max, tol, nls_method, cg_method);
% [data.xsol, data.logg, flag] = test_model.solve_nls;            
% 
% 
% data.post_dataprocessing(test_model.r);
% % iter Information processing
% %% results
% data.post_dataprocessing(test_model.r);
% 
% 
% %%
% nls_method  = 'preconditioned Levenberg-Marquardt';
% cg_method  = 'CG-Steihaug';
% test_model.options = nlsOption(iter_max, tol, nls_method, cg_method);
% [data.xsol, data.logg, flag] = test_model.solve_nls;                    
% 
% data.post_dataprocessing(test_model.r);
% 
% %%
% nls_method  = 'Gauss-Newton';
% cg_method  = 'preconditioned CG-Steihaug';
% test_model.options = nlsOption(iter_max, tol, nls_method, cg_method);
% [data.xsol, data.logg, flag] = test_model.solve_nls;                    
% track = [data.x0,data.logg.xk];
% plot(track(1,:),track(2,:))
% 
% data.post_dataprocessing(test_model.r);
% 
% 


close all
clear 
addpath(genpath('../data/'));
addpath(genpath('../class_def/'));
addpath(genpath('../nls_solvers/'));
%%
n =10;
cond_A = zeros(n,1);
e_inv = zeros(n,1);
t_inv = zeros(n,1);
e_pinv = zeros(n,1);
t_pinv = zeros(n,1);

e_cg = zeros(n,1);
t_cg = zeros(n,1);


i_cg = zeros(n,1);

e_pcg = zeros(n,1);
t_pcg = zeros(n,1);

i_pcg = zeros(n,1);
N     = zeros(n,1);
options.diagcomp = 10^-6;
% options.type = 'ict';
% options.droptol = 1e-2;
% options.michol   = 'on';
tol   = 10^-5;
for i = 1:n
    % create problem
    N(i)             = 2^i;
    [A,b,xref]       = ill_conditioned_example(N(i));
    
    cond_A(i)        = cond(A);
    % inv
    tic
    x                = inv(A)*b;
    t_inv(i)         = toc;
    e_inv(i)         = norm(xref - x,inf);
    % pinv
    
    tic
    x                = pinv(A,1e-5)*b;
    t_pinv(i)        = toc;
    e_pinv(i)        = norm(xref - x,inf);
    % cg
    tic
    [x_cg,i_cg(i)]   = pcg_steihaug(A,b,[],N(i),[]);
    t_cg(i)          = toc;
    e_cg(i)          = norm(xref - x_cg,inf);
    
    
    
    % pcg
    A                = sparse(A);
    tic
    C_T = ichol(A,options);
    [x_pcg,i_pcg(i)] = pcg_steihaug(A,b,[],N(i),C_T);    
    t_pcg(i)         = toc;
    e_pcg(i)         = norm(xref - x_pcg,inf);
    % pcg-2
%     dk                = diag(ichol(A));
%     dk(dk==0)        = 1e-5;
%     Dk               = diag(dk);
%     [x_cg,i_cg(i)] = pcg_steihaug(A,b,[],N(i),Dk);    
%     t_cg(i)         = toc;
%     e_cg(i)         = norm(xref - x_cg,inf);

% pcg-3
%     [x_cg,i_cg(i)]   =  pcg(A,b,tol,N(i),C_T);
%     t_cg(i)          = toc;
%     e_cg(i)          = norm(xref - x_cg,inf);
   
    
end
%%
plot_ill_condition_result(t_inv,t_pinv,t_pcg,t_cg, e_inv,e_pinv,e_cg,e_pcg,i_cg,i_pcg,N,cond_A)

function plot_ill_condition_result(t_inv,t_pinv,t_pcg,t_cg,e_inv,e_pinv,e_cg,e_pcg,i_cg,i_pcg,N,cond_A)

figure('Name','Condition Number for Inversion')
semilogy(cond_A)
xlabel('$\mathrm{Dimesion\;of\;Hilbert\;Martix}$','fontsize',12,'interpreter','latex')
ylabel('$||x^*-x^{\mathrm{sol}}||_\infty$','fontsize',12,'interpreter','latex')
set(gca, 'XTickLabel',N)
% xlim(x_bound)
grid on 



figure('Name','Error of Ill Conditioned Problem')
semilogy([e_pcg,e_cg,e_inv,e_pinv])
legend({'Preconditioned CG','CG','inv()','pinv()'},'fontsize',12,'interpreter','latex','Location','northwest')
xlabel('$\mathrm{Dimesion\;of\;Hilbert\;Martix}$','fontsize',12,'interpreter','latex')
ylabel('$||x^*-x^{\mathrm{sol}}||_\infty$','fontsize',12,'interpreter','latex')
set(gca, 'XTickLabel',N)
% xlim(x_bound)
ylim([10^-16,10^15])
grid on 

figure('Name','Computing time')
subplot(2,1,1)
semilogy([t_pcg,t_cg,t_inv,t_pinv])
legend({'Preconditioned CG','CG','inv()','pinv()'},'fontsize',12,'interpreter','latex','Location','northwest')
xlabel('$\mathrm{Dimesion\;of\;Hilbert\;Martix}$','fontsize',12,'interpreter','latex')
ylabel('$\mathrm{Computing\;Time[s]}$','fontsize',12,'interpreter','latex')
set(gca, 'XTickLabel',N)
grid on

subplot(2,1,2)
plot([i_pcg,i_cg])
xlabel('$\mathrm{Dimesion\;of\;Hilbert\;Martix}$','fontsize',12,'interpreter','latex')
ylabel('$\mathrm{Number\;of\;Iterations}$','fontsize',12,'interpreter','latex')
set(gca, 'XTickLabel',N)
grid on 

end


function [x,i] = pcg_steihaug(A,b,toltol,maxit,C_T)
% lin
% check alg options
if nargin < 3 || isempty(toltol), toltol = 1e-12; end
if nargin < 4 || isempty(maxit), maxit = min(20,length(b)); end
% if nargin <5 || isempty(delta), delta = 1; end
% intial
x = zeros(size(b));
r = -b;
if ~isempty(C_T)
    % preconditioning activated
    y = C_T'\(C_T\r);
    p = - y;
    rr = r'*y;    
else
    % without preconditioning
    p = -r;
    rr = r'*r;
end
Ap = A*p;
i = 1;
while (max(abs(r))>eps) &&  i<=maxit

    rho = rr/(p'*Ap);       % one-dim minimizer
    x   = x + rho*p;        % update state

    r   = r + rho*Ap;       % update residual
    if ~isempty(C_T)
        % preconditioning activated
        y   = C_T'\(C_T\r);
        rr_new = r'*y;
    else
        % without preconditioning
        rr_new = r'*r;
    end
    if rr_new<toltol
        % 3. terminate when reach CG solution within trust region
        return
    end
    beta   = rr_new/rr;      % update the parameter to ensure conjugate 
    rr  = rr_new;
    if ~isempty(C_T)
        % preconditioning activated
        p   = -y + beta*p;
    else
        % without preconditioning
        p   = -r + beta*p;
    end    
    Ap  = A*p;
    i   = i+1;
end


end


function [A,b,xref] = ill_conditioned_example(n)
    % create a symmetric ill-conditioned matrix H
    A    = hilb(n);
    A(:,end) = 0;
    A(end,:) = 0;
    A(end,end) = 1e-10;
%     A    = A'*A;
    sigma = 3;
    random_vector = (2*rand(n,1)-1)*10;
    xref=random('Normal',random_vector,sigma*ones(n,1));  
    b    = A*xref;
end
