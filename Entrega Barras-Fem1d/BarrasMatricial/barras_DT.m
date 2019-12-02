function [Def,Ten] = barras_DT(xnode,icone,propiedades, model,U)
    
    for barra = 1 : model.nbarras
        indicesNodos = icone(barra,:);
    
        nodo1 = xnode(indicesNodos(1), :);
        nodo2 = xnode(indicesNodos(2), :);
        L = norm(nodo2 - nodo1);
        
        #a continuacion se obtienen los indices globales de las deformaciones de los nodos
        indx = [];
        for nodo = indicesNodos
            indx = [indx nodo*2-1 nodo*2];
        end
    
        alpha = atan2(nodo2(2) - nodo1(2), nodo2(1) - nodo1(1));
        c = cos(alpha);
        s = sin(alpha);
    
        C = 1/L * [ -c -s c s];
        
        deformacion = C * U(indx);
        Def(barra) = deformacion;
        
        E = propiedades(barra).E;
        Ten(barra) = E * deformacion;
    end
end
