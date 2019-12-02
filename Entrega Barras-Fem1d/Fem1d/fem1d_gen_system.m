function K = fem1d_gen_system(K,xnode,icone,model)
    % Matrices and vectors elementary -- assembly
    for e=1:model.nelem
        
        elem = icone(e, :);
        
        % Saco las filas que corresponde a los nodos que me interesan
        % elem contiene los indices de los nodos
        nodes = xnode(elem);
        
        % Local Matrices and Vectors
        localK = fem1d_genK(nodes,model.young,model.area);
        
        % Assembly
        % Ver a mano como funciona esto
        indx = [elem(1) elem(2)];
        
        K(indx,indx) += localK;
    end
end

