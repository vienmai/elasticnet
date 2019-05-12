function [out] =  ProxGD_Solver(x0, f_eval, g_eval, Hess_opt, boot_method, opts)

if(~isempty(boot_method))
    boot_method(x0,f_eval, g_eval, Hess_opt,opts);
end

x = x0;
f_full = @(x) f_eval(x, 1:opts.numdata);
f0 = f_full(x0);
f = f0;
iteration =1;
cur_time =0; xdiff =0;
outer_count = 1;

% % Information for plotting
if(opts.plotting)
    errors = f0;%1;
    times = 0;
    DATA.datapasses = 0;
    DATA.datapasses_products =0;
end

while (iteration <= opts.max_iterations)
    tic; % START timing
    
    DATA.grad  = g_eval(x,1:opts.numdata);
    
    x = prox_l1(x - opts.step_parameter*DATA.grad, opts.lambda*opts.step_parameter) ;
    
    cur_time =cur_time+ toc;  % END Timing iteration
    
    
    f =  f_full(x)
    %         opts.x_nnz(outer_count) = nnz(x);
    %         nnz(x)
    % Information for plotting
    
    if(isnan(f) || isinf(f))
        out.stopping_flag = 'NaN';
        break;
    end
    if(opts.plotting)
        errors(end+1) = f;% f/f0;
        times(end+1) = cur_time;
    end
    iteration = iteration+1;
end % End of solver
%% Wrap up output
% Time the finish, remainder is for testing
out.cput = toc;
if (cur_time > opts.Timeout)
    out.stopping_flag = 'timeout';
elseif (isnan(f) || isnan(xdiff) || isinf(f) || isinf(xdiff))
    out.stopping_flag = 'NaN';
else
    out.stopping_flag = 'conv.';
end

out.name = DATA.name;
out.iteration = iteration;
out.x = x;
out.f = f;

if(opts.plotting)
    out.times =times;
    out.errors =errors;
    out.datapasses =DATA.datapasses;
end
%% Test the output
out.accuracy=  accurary_prediction( opts.X, opts.y, x);
%   out.val_accuracy = test_log_model(opts, out);
opts.X = [];
opts.y = [];
opts.x0 = [];
out.opts = opts;

%% Print conclusion
if opts.prnt
    iteration  =0;
    if(strcmp(opts.grad_type,'SGD'))
        fprintf('\n(beta,alpha) = (%f, %f)\n', opts.beta,opts.alpha );
    end
    print_iteration_info;
end
fprintf('%s total time %12.6f s', out.name,cur_time); display('   ');
end


