function g = ls_grad_sub(A, b, x, l2_reg, S)

g = ls_grad(A(S,:), b(S), x, l2_reg);

end