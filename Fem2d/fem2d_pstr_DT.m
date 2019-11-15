function [Def_prom,Ten_prom,Ten_VM] = fem2d_pstr_DT(xnode,icone,model,D,U)
    
    esTriangular = @(elemento) icone(elemento,4) == -1;
    
    Def = zeros(model.nnodes,4);
    
    for elemento = 1 : model.nelem
        
        if (esTriangular(elemento))
            indicesNodos = icone(elemento,1:3);
        else
            indicesNodos = icone(elemento,:);
        end
        
        indicesU = [];
        for indiceNodo = indicesNodos
            indicesU = [indicesU indiceNodo*2-1 indiceNodo*2];
        end
        
        nodos = xnode(indicesNodos,:);
        desplazamientos = U(indicesU);
        
        if (esTriangular(elemento))
            x = nodos(:, 1);
            y = nodos(:, 2);
            
            J = [x(2)-x(1)  y(2)-y(1);
                 x(3)-x(1)  y(3)-y(1)];
            DN = [-1 1 0; -1 0 1];
            V = inv(J)*DN;
            B = [V(1,1)     0     V(1,2)     0     V(1,3)     0   ;
                   0      V(2,1)    0      V(2,2)    0      V(2,3);
                 V(2,1)   V(1,1)  V(2,2)   V(1,2)  V(2,3)   V(1,3)];
            Def_xy = (B*desplazamientos)';
            
            for k = 1 : 3
                Def(indicesNodos(k),1:3) += Def_xy;
                Def(indicesNodos(k),4) += 1;
            end
        else 
            % cuatro puntos de Gauss con peso w=1
            p = [-1 -1;
                  1 -1;
                  1  1;
                 -1  1];
            
            for ii = 1 : 4
                    s = p(ii,1);
                    t = p(ii,2);
                    DNnum = [(-1+t)/4   (1-t)/4    (1+t)/4   (-1-t)/4;
                             (-1+s)/4   (-1-s)/4   (1+s)/4   (1-s)/4];
                    J = DNnum*nodos;
                    V = inv(J)*DNnum;
                    B = [V(1,1)     0     V(1,2)     0     V(1,3)     0     V(1,4)     0;
                           0      V(2,1)    0      V(2,2)    0      V(2,3)    0      V(2,4);
                         V(2,1)   V(1,1)  V(2,2)   V(1,2)  V(2,3)   V(1,3)  V(2,4)   V(1,4)];
                    
                    Def_xy = (B*desplazamientos)';
                    Def(indicesNodos(ii),1:3) += Def_xy;
                    Def(indicesNodos(ii),4) += 1;
            end
        end
    end
    
    Def_prom(:,1) = Def(:,1) ./ Def(:,4);
    Def_prom(:,2) = Def(:,2) ./ Def(:,4);
    Def_prom(:,3) = Def(:,3) ./ Def(:,4);
    
    for i = 1 : model.nnodes
        Ten_prom(i,1:3) = (D*Def_prom(i,1:3)')';
    end
    
    Ten_VM = Ten_prom(:,1).^2 - Ten_prom(:,1).*Ten_prom(:,2) + Ten_prom(:,2).^2 ...
            + 3*(Ten_prom(:,3).^2);
end 