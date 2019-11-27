function [F] = fem1d_pointload(F,pointload)
    
    for ii = 1 : size(pointload, 1)
        
        nodo = pointload(ii, 1);
        nodo = double(nodo);
        X = pointload(ii, 2);
        
        F(nodo) += X;
    end
end