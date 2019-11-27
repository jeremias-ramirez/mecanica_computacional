function [U,reaction,Def,Ten] = barras(xnode,icone,propiedades,Fixnodes,Pointload,model)
    % Initialize system matrices and vectors
    [K,F] = barras_initialize(model.nnodes);
    
    % Assembly system's matrix K and vector F
    [K,F] = barras_gen_system(K,F,xnode,icone,propiedades,model);
    
    % Point loads
    [F] = barras_pointload(F,Pointload);
    
    % Displacements Condition
    [Kd,Fd] = barras_fixnodes(K,F,Fixnodes);

    % Solve linear equations system
    U = Kd\Fd;
    
    % Compute stress and strains from displacements U
    [Def,Ten] = barras_DT(xnode,icone,propiedades,model,U);
    
    % Compute the reactions on the fixed nodes as a R = StifMat * u - F
    reaction = K*U;
end
