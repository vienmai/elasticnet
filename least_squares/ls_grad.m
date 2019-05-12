function g = ls_grad(A, b, x, l2_reg)
n=size(A,1); 
g=A'*(A*x-b)/n+l2_reg*x;   
end