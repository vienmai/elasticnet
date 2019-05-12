function Hv = ls_hessv(hess, v)
% Compute H*v

tmp = hess.V'*v;
Hv = hess.V*(hess.S*tmp) + hess.S(end,end)*(v-hess.V*tmp);

end