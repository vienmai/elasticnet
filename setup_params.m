function [opts, g_stoc_eval] = setup_params(opts,hess,method)

switch method
    case 'acc_scaled_prox_svrg'
        opts.num_epoch_per_outer_loop = 2;
        opts.cond='avg';
        [opts.L, opts.Lvec]= Lipschitz_eval(opts, hess);
        opts.mu = opts.l2_reg;
        opts.cvxcomb = sqrt(opts.mu/opts.L);
%         opts.S = min(floor(sqrt(opts.L/opts.mu)), floor(sqrt(opts.numdata)));
%         opts.S = floor(sqrt(opts.numdata));
        opts.S = 500;
        opts.inner_iters = floor(opts.num_epoch_per_outer_loop*(opts.numdata/opts.S));
        fprintf('*********************\n')
        fprintf('\n')
        fprintf('Lipschitz: %f \t Batch-size: %6d \n', opts.L, opts.S);
        fprintf('\n')
    case 'katyusha'
        opts.num_epoch_per_outer_loop = 2;
        opts.cond='l2_avg';
        [opts.L, opts.Lvec]= Lipschitz_eval(opts, hess);
        opts.mu = opts.l2_reg;
                
        opts.S = floor(sqrt(opts.numdata));
        opts.inner_iters = floor(opts.num_epoch_per_outer_loop*(opts.numdata/opts.S));
        
        [opts]=katyusha_params(opts);

        fprintf('*********************\n')
        fprintf('\n')
        fprintf('Lipschitz: %f \t Batch-size: %6d \n', opts.L, opts.S);
        fprintf('\n')
    case 'prox_svrg'
        opts.num_epoch_per_outer_loop = 2;
        opts.cond='l2_avg';
        [opts.L, opts.Lvec]= Lipschitz_eval(opts, hess);
        opts.mu = opts.l2_reg;
        
        opts.S = floor(sqrt(opts.numdata));
        opts.inner_iters = floor(opts.num_epoch_per_outer_loop*(opts.numdata/opts.S));
        
        fprintf('*********************\n')
        fprintf('\n')
        fprintf('Lipschitz: %f \t Batch-size: %6d \n', opts.L, opts.S);
        fprintf('\n')
    otherwise
        error(['Unknown method ' method]);
end

g_stoc_eval = @(x,S)ls_grad_sub_non_unif(opts.A, opts.b, x, opts.l2_reg, opts.L, opts.Lvec, S); % nonuniform sampling 

end