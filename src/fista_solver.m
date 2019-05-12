function [out] =  fista_solver(x0, L, opts, max_iter)

%Params for FISTA
f = @(x)0.5*(norm(opts.A*x-opts.b)^2)/opts.numdata + 0.5*opts.l2_reg*norm(x)^2;
g = @(x)opts.A'*(opts.A*x-opts.b)/opts.numdata + opts.l2_reg*x;
h = @(x)norm(x,1);
prox_h = @(x,a)prox_l1(x,a);

par.max_iter = max_iter;
par.Lstart = L;
par.const_flag = true;
par.print_flag = true;
par.eps = 1e-15;

[x, fmin, parout_fista] = fista(f, g, h, prox_h, opts.l1_reg, x0, par);

%% Wrap up output
out.x = x;
out.f = fmin;
out.times =parout_fista.times;
out.errors =parout_fista.funValVec;
out.passes =parout_fista.iterNum;
end


