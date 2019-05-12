function [out] =  acc_scaled_prox_svrg(x0, g_eval, g_stoc_eval, hess, opts)

x = x0;
z = x;
w0 = x0;

f_full = @(x) ls_eval(opts.A, opts.b, x, opts.l2_reg, opts.l1_reg);
f0 = f_full(x0);

iter = 0;
cur_time = 0;
outer_count = 0;

% % Information for plotting
errors = f0;
times = 0;
DATA.passes = 0;

if(opts.is_print)
    fprintf('*********************\n');
    fprintf('HSke\n') ;
    fprintf('*********************\n')
    fprintf('#iter       \tfun. val.\n');
end

while (iter <= opts.max_iter*opts.inner_iters)
        
    % initilize the subproblem solver
    if iter ~= 0
        w0 = warm_start(x, y_, opts.step_size/hess.S(1,1), hess, opts);
    end
%     w0 = x;
    tic; % START timing
    
    % main iterates
    y = x/(1+opts.cvxcomb) + opts.cvxcomb*z/(1+opts.cvxcomb);
    step_scaled_prox_svrg
    Hinv = ls_hessinvv(hess, stoc_grad);
    y_ = y - opts.step_size*Hinv;
    x = scaled_l1_prox(w0, y_, hess, opts, opts.solver);    
    z = z + opts.cvxcomb*(y-z) - (opts.cvxcomb/opts.mu/opts.step_size)*(y-x);
    
    cur_time = cur_time + toc;  % END Timing iteration
    
    if(mod(iter,opts.inner_iters)==0)
        f = f_full(x);
        outer_count = outer_count + 1;
        if(isnan(f) || isinf(f)||f>1e2)
            out.stopping_flag = 'NaN';
            fprintf('I am breaking the loop.\n')
            break;
        elseif f>errors(end)
            out.stopping_flag = 'noDesc';
            break;    
        end
        errors(end+1) = f;
        times(end+1) = cur_time;
        if (opts.is_print)
            fprintf('%6d \t    %12f \t',outer_count, f);
            fprintf('\n') ;
        end
    end
    iter = iter+1;
    if(cur_time > opts.Timeout)
        fprintf('Time up')
        break;
    end
    if( DATA.passes(end) > opts.totalpasses)
        break;
    end
end % End of solver

if (opts.is_print)
    fprintf('----------------------------------\n') ;
    fprintf('Optimal value = %15f \n',f) ;
    fprintf('----------------------------------\n') ;
else
    fprintf('Optimal value = %15f \n',f) ;
end

%% Wrap up output
out.cput = toc;
if (cur_time > opts.Timeout)
    out.stopping_flag = 'timeout';
elseif (isnan(f) || isinf(f) || f > 1e2)
    out.stopping_flag = 'NaN';
else
    out.stopping_flag = 'OK';
end

out.iter = iter;
out.x = x;
out.f = f;
out.times =times;
out.errors =errors;
out.passes =DATA.passes;
end


