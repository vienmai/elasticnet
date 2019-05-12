function [OUT, beststep] =  optimal_step_size(opts, hess, g_eval, g_stoc_eval, grid, method)
OUT={};
minf = inf;
beststep = grid(end);
opts.get_optimal_step_size =0;
fprintf('*********************\n')
for gg = 1: length(grid)
    try
        switch method
            case 'acc_scaled_prox_svrg'
                fprintf('\n');
                display(['Trying stepsize: ' num2str(grid(gg))]);
                opts.step_size = grid(gg);
                [out] = acc_scaled_prox_svrg(opts.x0, g_eval, g_stoc_eval, hess, opts);
            case 'katyusha'
                fprintf('\n');
                display(['Trying stepsize: ' num2str(grid(gg))]);
                opts.katy_step_size = grid(gg);
                [out] =  katyusha1(opts.x0, g_eval, g_stoc_eval, opts);
            case 'prox_svrg'
                fprintf('\n');
                display(['Trying stepsize: ' num2str(grid(gg))]);
                opts.step_size = grid(gg);
                [out] = prox_svrg(opts.x0, g_eval, g_stoc_eval, opts);
            otherwise
                error(['Unknown method ' method]);
        end
    catch something_messed_up
        continue;
    end
    if(strcmp(out.stopping_flag,'NaN') || (max(out.errors) > out.errors(1)))
        continue;
    elseif (out.errors(end) < minf)        
        minf = out.errors(end);
        beststep = grid(gg);
        outbest = out;
    end
end
OUT=outbest;
display(['Best step: ' num2str(beststep)]);
end


