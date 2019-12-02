function [K,F] = barras_gen_system(K,F,xnode,icone,propiedades,model)
    
    for barra = 1 : model.nbarras
        indicesNodos = icone(barra,:);
        
        % Saco las filas que corresponde a los nodos que me interesan
        % elem contiene los indices de los nodos
        nodes = xnode(indicesNodos,:);
        
        % Local Matrices and Vectors
        localK = barras_genK(nodes,propiedades(barra));
        
        % Assembly
        indx = [];
        for nodo = indicesNodos
            indx = [indx nodo*2-1 nodo*2];
        end
        
        K(indx,indx) += localK;
    end
    
end
