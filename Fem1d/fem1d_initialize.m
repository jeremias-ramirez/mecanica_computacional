function [K,F] = fem1d_initialize(nnodes)
    % Initialize system matrices and vectors
    dof = nnodes;   % degree of freedom
    K = zeros(dof,dof);
    F = zeros(dof,1);
end

