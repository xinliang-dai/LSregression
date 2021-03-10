function x0 = vary_initial_point(x_sol,sigma)
    
    Nx            =  numel(x_sol);
%     random_vector =  (2*rand(Nx,1)-1)*radiu;
    x0            =  random('Normal',x_sol,sigma*ones(Nx,1));  
end