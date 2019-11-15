function [K,F] = fem2d_pstr_gen_system(K,F,D,xnode,icone,model)
    % Matrices and vectors elementary -- assembly
    for e=1:model.nelem
        if (icone(e,4) == -1)
            elem = icone(e,1:3);
        else
            elem = icone(e,:);
        end
        
        % Saco las filas que corresponde a los nodos que me interesan
        % elem contiene los indices de los nodos
        nodes = xnode(elem,:);
        
        % Local Matrices and Vectors
        localK = fem2d_pstr_genK(nodes,D,model.thick);
        localF = fem2d_pstr_genF(nodes,model.gravity,model.thick);
        
        % Assembly
        % Ver a mano como funciona esto
        indx = [];
        for i = 1 : length(elem)
            indx = [indx elem(i)*2-1 elem(i)*2];
        end
        
        K(indx,indx) = K(indx,indx) + localK;
        F(indx) = F(indx) + localF;
    end
end

