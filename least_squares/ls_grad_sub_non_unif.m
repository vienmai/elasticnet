function g = ls_grad_sub_non_unif(A, b, x, l2_reg, L, Lvec, S)

g = ls_grad_non_unif(A(S,:), b(S), x, l2_reg, L, Lvec(S));

end