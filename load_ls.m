function [opts, g_eval, g_stoc_eval] =  load_ls(datapath, tol, opts)

if(strcmp(opts.dataset,'cina0'))
    load('./data/cina0.mat');
    opts.A = A;
    opts.b = b;
else
    datafile = [datapath '/' opts.dataset];
    [opts.b, opts.A] = libsvmread(datafile);
end
sX = size(opts.A);
opts.numfea = sX(2);
opts.numdata = sX(1);
opts.x0  = zeros(opts.numfea,1);
opts.tol = tol;
opts.Timeout = inf;
opts.max_iter = inf;
opts.l2_reg = 1e-3;
opts.l1_reg = 1e-3;

g_eval = @(x,S)ls_grad_sub(opts.A, opts.b, x, opts.l2_reg, S);  %used for fullgrad and uniform sampling
end
