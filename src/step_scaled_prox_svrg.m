% Outer loop step
if(mod(iter,opts.inner_iters)==0)
    inner_count =0;
    % sampling with replacement
    rand_samples = randsample(opts.numdata, opts.num_epoch_per_outer_loop*opts.numdata, true, opts.Lvec);
%     rand_samples = randsample(opts.numdata, opts.inner_iters*opts.S, true);
    xref=x;
    y = xref;
    z = xref;
    % Calculate full gradient
    full_grad =  g_eval(xref,1:opts.numdata);
    % Calculate the data passes
    DATA.passes(end+1) = DATA.passes(end) + 1;
    if(iter ~= 0)
        DATA.passes(end) = DATA.passes(end) + opts.inner_iters*opts.S/opts.numdata;
    end
end
% Inner loop step

if (inner_count~=opts.inner_iters-1)
    grad_sample = rand_samples(inner_count*opts.S+1:(inner_count+1)*opts.S);
else
    grad_sample = rand_samples(inner_count*opts.S+1:end);
end
stoc_grad =  g_stoc_eval(y, grad_sample) - g_stoc_eval(xref,grad_sample) + full_grad;
inner_count = inner_count + 1;


