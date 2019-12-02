function [U,reaction] = fem1d(xnode,icone,Fixnodes,Sideload,Pointload,model)
    % Initialize system matrices and vectors
    [K,F] = fem1d_initialize(model.nnodes);
    
    % Assembly system's matrix K and vector F
    K = fem1d_gen_system(K,xnode,icone,model);
    
    % Surface tractions
    F = fem1d_sideload(F,Sideload,xnode);
    
    % Point loads
    F = fem1d_pointload(F,Pointload);
    
    % Displacements Condition
    [Kd,Fd] = fem1d_fixnodes(K,F,Fixnodes);

    % Solve linear equations system
    U = Kd\Fd;
    
    % Compute the reactions on the fixed nodes as a R = StifMat * u - F
    reaction = K * U - F;
end
