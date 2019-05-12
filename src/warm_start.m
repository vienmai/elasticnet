function z = warm_start(x, u, ws_stepsize, hess, opts)

tmp = ls_hessv(hess, x-u)/opts.step_size;

z = prox_l1(x - ws_stepsize*tmp, ws_stepsize*opts.l1_reg);