function [out] =  bcd(x0, opts)
 
x = x0;

f_full = @(x) ls_eval(opts.A, opts.b, x, opts.l2_reg, opts.l1_reg);
f0 = f_full(x0);

r = opts.b - (opts.A)*x;

cur_time = 0;
outer_count = 0;

% % Information for plotting
errors = f0;
times = 0;

if(opts.is_print)
    fprintf('*********************\n');
    fprintf('BCD\n') ;
    fprintf('*********************\n')
    fprintf('#iter       \tfun. val.\n');
end

while (outer_count <= opts.max_iter)
    tic; % START timing
    
    for idx=1:length(x)
        
        normAi = norm(opts.A(:,idx));
        beta = normAi^2/opts.numdata + opts.l2_reg;
        
        Aixi = opts.A(:,idx)*x(idx);
        
        bcd_step = (opts.A(:,idx)'*r/opts.numdata + x(idx)*(normAi^2)/opts.numdata)/beta;
        
        x(idx) = prox_l1(bcd_step, opts.l1_reg/beta); 
        
        r = r + Aixi - opts.A(:,idx)*x(idx);
        
    end
    
    cur_time = cur_time + toc;  % END Timing iteration

    f = f_full(x);
    outer_count = outer_count + 1;
    if(isnan(f) || isinf(f) || f>1e2)
        out.stopping_flag = 'NaN';
        break;
    end
    errors(end+1) = f;
    times(end+1) = cur_time;
    if (opts.is_print)
        fprintf('%6d \t    %12f \t',outer_count, f);
        fprintf('\n') ;
    end

    if(cur_time > opts.Timeout)
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

out.iter = outer_count;
out.x = x;
out.f = f;
out.times =times;
out.errors =errors;
end


