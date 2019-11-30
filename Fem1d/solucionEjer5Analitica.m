function y = solucionEjer5Analitica(x, Q, L, E, A)
    
    y = zeros(length(x), 1);
    
    if (length(findsym(x)) != 0)
        y = sym(y);
    end

    indiceL = find(x - L < 1e-5);
    indiceL = double(indiceL);
    
    if isempty(indiceL)
        error("No se encontró elemento con valor L");
    end
    
    indicesMenor = 1 : indiceL - 1;
    y(indicesMenor) = 2/9 * Q/(E*A) * x(indicesMenor);
    
    indicesMayor = indiceL : length(x);
    xM = x(indicesMayor);
    y(indicesMayor) = 1/36 * Q*L/(E*A) * (3 - xM/L + 9 * (xM/L).^2 - 3 * (xM/L).^3);
end
