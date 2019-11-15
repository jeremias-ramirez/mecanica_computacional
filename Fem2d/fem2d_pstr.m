function [U,reaction,Def,Ten,Ten_VM] = fem2d_pstr(xnode,icone,Fixnodes,Sideload,Pointload,model)
    % Initialize system matrices and vectors
    [K,F] = fem2d_pstr_initialize(model.nnodes);
    
    % Obtain constitutive matrix
    [D] = fem2d_pstr_const_mat(model);
    
    % Assembly system's matrix K and vector F
    [K,F] = fem2d_pstr_gen_system(K,F,D,xnode,icone,model);
    
    % Surface tractions
    [F] = fem2d_pstr_sideload(F,Sideload,xnode,model.thick);
    
    % Point loads
    [F] = fem2d_pstr_pointload(F,Pointload,xnode,icone);
    
    % Displacements Condition
    [Kd,Fd] = fem2d_pstr_fixnodes(K,F,Fixnodes);

    % Solve linear equations system
    U = Kd\Fd;
    
    % Compute stress and strains from displacements U
    [Def,Ten,Ten_VM] = fem2d_pstr_DT(xnode,icone,model,D,U);
    
    % Compute the reactions on the fixed nodes as a R = StifMat * u - F
    [reaction] = fem2d_pstr_reaction(K,F,U,Fixnodes);
end
