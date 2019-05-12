function x = scaled_l1_prox(z0, u, hess, opts, method)

subcond = hess.S(1,1)/hess.S(end,end);

switch method
    case 'fista'

        f = @(x)0.5*(x-u)'*(ls_hessv(hess,x-u))/opts.step_size;
        g = @(x)ls_hessv(hess,x-u)/opts.step_size;
        h = @(x)norm(x,1);
        prox_h = @(x,a)prox_l1(x,a);
        
        par.max_iter = sqrt(subcond)*log10(subcond);
        par.Lstart = hess.S(1,1)/opts.step_size;
        par.print_flag = false;
        par.eco_flag = true;
        par.const_flag = true;
        par.eps = 1e-15;

        x = fista_(f, g, h, prox_h, opts.l1_reg, z0, par);
       
    case 'bcd'

        max_iter = 0; 
        x =  bcd_ls(z0, hess, u, opts.step_size*opts.l1_reg, max_iter);        
    
    otherwise
        error(['Unknown method ' method]);
end

