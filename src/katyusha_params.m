function [opts]=katyusha_params(opts)

if opts.S==1
    opts.katy_tau_1 = min(0.5, sqrt(opts.inner_iters*opts.mu/3/opts.L));
    opts.katy_tau_2 = 0.5;
    opts.katy_alpha = 1/3/opts.katy_tau_1/opts.L;
else
    opts.katy_tau_2 = min(0.5*opts.L/opts.bestL/opts.S,0.5);
%     opts.katy_tau_2 = 0.5*opts.S;
    
    if opts.bestL <= opts.L*opts.inner_iters/opts.S
        opts.katy_tau_1 = min(sqrt(8*opts.S*opts.inner_iters*opts.mu/3/opts.L)*opts.katy_tau_2,opts.katy_tau_2);
        estL = 0.5*opts.L/opts.S/opts.katy_tau_2;
    else
        opts.katy_tau_1 =  min(sqrt(2*opts.mu/3/opts.bestL),0.5/opts.inner_iters);
        estL=opts.bestL;
    end
    opts.katy_alpha = 1/3/opts.katy_tau_1/estL;
end