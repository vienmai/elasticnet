clc; 
clear;

opts.dataset = 'gisette_scale';
datapath = './data';
tol = 0;
[opts, g_eval] =  load_ls(datapath, tol, opts);

opts.totalpasses = 100;
opts.Timeout = 400;

%% Hessian approximation
opts.approx_rank = 40;
oversamp = min(5,abs(opts.numfea-opts.approx_rank));
[hess.U, hess.S, hess.V] = bksvd(opts.A/sqrt(opts.numdata), opts.approx_rank, 2, opts.approx_rank+oversamp);
hess.S = (hess.S)^2 + opts.l2_reg*eye(size(hess.S));
hess.sqrtS = sqrtm(hess.S);
opts.bestL = hess.S(1,1);

%% Grid of steps

% grid = [1e-2, 2e-2, 5e-2, 1e-1, 2e-1, 5e-1, 1e0, 2e0, 5e0, 1e1, 2e1, 5e1, 1e2, 2e2, 5e2];

%% HSKE

method='acc_scaled_prox_svrg';
opts.is_print = 1;

opts.solver = 'fista';

[opts, g_stoc_eval] = setup_params(opts,hess,method);
[out_hse, opts.step_size] =  optimal_step_size(opts, hess, g_eval, g_stoc_eval, grid/opts.L, method);


% opts.step_size = 0.674483721748090;
  
% [out_hse] = acc_scaled_prox_svrg(opts.x0, g_eval, g_stoc_eval, hess, opts);

%% BCD

opts.is_print = 1;
opts.max_iter = 1000;
[out_bcd] = bcd(opts.x0, opts);


%% FISTA

max_iter = 1500;
[out_fista] =  fista_solver(opts.x0, hess.S(1,1), opts, max_iter);

%% Katyusha

method='katyusha'; % acc_scaled_prox_svrg  katyusha   prox_svrg;
opts.is_print = 0;

[opts, g_stoc_eval] = setup_params(opts,hess,method);
[out_, katy_best_step] = optimal_step_size(opts, hess, g_eval, g_stoc_eval, grid*opts.katy_alpha, method);
opts.katy_step_size = katy_best_step; 
[out_] = katyusha1(opts.x0, g_eval, g_stoc_eval, opts);
% 
%% Prox-SVRG
method='prox_svrg';
opts.is_print=0;

[opts, g_stoc_eval] = setup_params(opts,hess,method);
[out_svrg, svrg_best_step] =  optimal_step_size(opts, hess, g_eval, g_stoc_eval, grid/opts.L, method);
opts.step_size = svrg_best_step; %0.05
[out_svrg] = prox_svrg(opts.x0, g_eval, g_stoc_eval, opts);

%% Plotting

optval= min([out_hse.errors(end), out_fista.errors(end), out_bcd.errors(end)]);

semilogy(out_hse.passes, out_hse.errors-optval,'-b')
hold  on

semilogy([0 1:out_fista.passes], out_fista.errors-optval,'c')
hold on

semilogy(out_katy.passes, out_katy.errors-optval,'-r')
hold  on

semilogy(out_svrg.passes, out_svrg.errors-optval,'k')

semilogy([0 1:out_bcd.iter], out_bcd.errors-optval,'-g')
 

xlabel('runtime [s]')
ylabel('residual')
title('australian')
legend('Hske','FISTA','Katyusha','Prox-SVRG', 'BCD')

 
%% Runtime plots
figure(2)

semilogy(out_hse.times, out_hse.errors-optval,'-b')
hold  on

semilogy(out_bcd.times, out_bcd.errors-optval,'-g')
hold  on


