function g = ls_grad_non_unif(A, b, x, l2_reg, L, Lvec)

n=size(A,1); 
tmp=A'*(diag(Lvec)\(A*x-b))+l2_reg*x*sum(1./Lvec);   
g = L*tmp/n;

end

