function [out] =  bcd_ls(x0, hess, u, l1_reg, max_iter)
x = x0;
v = u-x0;
tmp = hess.V'*v;
c = sqrt(hess.S(end,end));
% r = hess.V*(hess.sqrtS*tmp) + c*(v-hess.V*tmp); % r0 = b-Ax_0 = A(u-x_0); A= sqrt(H)
r = hess.V*(hess.sqrtS*tmp - c*tmp) + c*v;
E = eye(length(x));

outer_count = 0;
while (outer_count <= max_iter)   
    for idx=1:length(x)        
        Ai = hess.V*(hess.sqrtS*hess.V(idx,:)'- c*hess.V(idx,:)') + c*E(:,idx); 
        normAi = norm(Ai);
        Aixi = Ai*x(idx);
        bcd_step = Ai'*r/normAi^2 + x(idx);
        x(idx) = prox_l1(bcd_step, l1_reg/normAi^2); 
        r = r + Aixi - Ai*x(idx);
    end
    outer_count = outer_count + 1;
end % End of solver
out = x;
end
