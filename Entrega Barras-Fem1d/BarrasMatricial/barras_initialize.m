function [K,F] = barras_initialize(nnodes)
    % Initialize system matrices and vectors
    dof = nnodes*2;   % degree of freedom
    K = zeros(dof,dof);
    F = zeros(dof,1);
end