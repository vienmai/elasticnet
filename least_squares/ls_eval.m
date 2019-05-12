function out = ls_eval(A, b, x, l2_reg, l1_reg)
 n=size(A,1); 
 out =  (0.5*norm(A*x-b)^2)/n + 0.5*l2_reg*norm(x)^2 + l1_reg*norm(x,1);
end