function [F] = fem2d_pstr_pointload(F,pointload,xnode,icone)
    
    esTriangular = @(elemento) icone(elemento,4) == -1;
    
    for n = 1 : size(pointload,1)
        elemento = pointload(n,1);
        Fx = pointload(n,2);
        Fy = pointload(n,3);
        xp = pointload(n,4);
        yp = pointload(n,5);
        
        if esTriangular(elemento)
            indicesNodos = icone(elemento,1:3);
        else
            indicesNodos = icone(elemento,:);
        end
        
        nodos = xnode(indicesNodos,:);
        N = fem2d_pstr_blerp(nodos,xp,yp);
        indx = [];
        f = [];
        
        for i = 1 : length(indicesNodos)
            indx = [indx indicesNodos(i)*2-1 indicesNodos(i)*2];
            f = [f N(i)*Fx N(i)*Fy];
        end
        
        F(indx) += f';
    end
end