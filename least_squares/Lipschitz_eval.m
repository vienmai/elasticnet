function [L, Lvec]=Lipschitz_eval(opts, hess) 

switch opts.cond
    case 'max' 
        for i=1:opts.numdata
            tmp = ls_hessinvv(hess,opts.A(i,:)');
            Lvec(i)= opts.A(i,:)*tmp;
        end
        L = max(Lvec);
    case 'avg' 
        for i=1:opts.numdata
            tmp = ls_hessinvv(hess,opts.A(i,:)');
            Lvec(i)= opts.A(i,:)*tmp;
        end
        L = sum(Lvec)/opts.numdata;
    case 'l2_max'
        Lvec = full(sum((opts.A).^2,2));
        L = max(Lvec);
    case 'l2_avg'
        Lvec = full(sum((opts.A).^2,2));
        L = sum(Lvec)/opts.numdata;
    otherwise
        display('Choose condition number type max or avg');
        error(['Unknown cond type ' opts.cond]);
end
