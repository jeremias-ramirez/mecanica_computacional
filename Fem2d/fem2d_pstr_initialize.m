function [K,F] = fem2d_pstr_initialize(nnodes)
    % Initialize system matrices and vectors
    dof = nnodes*2;   % degree of freedom
    K = sparse(dof,dof);
    F = sparse(dof,1);
end

