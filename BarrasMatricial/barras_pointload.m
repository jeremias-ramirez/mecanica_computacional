function [F] = barras_pointload(F,Pointload)
    
    for ii = 1 : size(Pointload, 1)
        nodo = Pointload(ii, 1);
        indexU = nodo*2 - 1;
        indexV = nodo*2;
        
        Fx = Pointload(ii, 2);
        Fy = Pointload(ii, 3);
        
        F(indexU) += Fx;
        F(indexV) += Fy;
    end
end
