function [F] = fem1d_sideload(F,Sideload,xnode)
% Descripción: módulo para calcular y ensamblar las contribuciones de pares de nodos
% (lados) donde se aplican cargas distribuidas.

% Las cargas distribuidas tienen forma b = alpha + beta * x;
% Entrada:
% * F: vector de fuerzas.
% * sideload: matriz con la información sobre elementos con cargas distribuidas. 
%   - Columnas 1-2: nodos del elemento en el cual actúa la carga.
%   - Columna 3: alpha.
%   - Columna 4: beta.
% * xnode: matriz de nodos con pares (x,y) representando las coordenadas de cada nodo
%   de la malla.

% Salida:
% * F: vector de fuerzas. Presenta modificaciones luego de aplicar la condición de borde.
% ----------------------------------------------------------------------

    for i = 1 : size(Sideload, 1)
        elem = Sideload(i, 1:2);
        alpha = Sideload(i, 3);
        beta = Sideload(i, 4);
        x1 = xnode(elem(1));
        x2 = xnode(elem(2));
        L = x2 - x1;
        
        F1 = 1/L * (alpha*(L^2)/2 + beta * (x1^3/3 - 0.5 * x1^2 * x2 + 1/6 * x2^3));
        F2 = 1/L * (alpha*(L^2)/2 + beta * (x1^3/6 - 0.5 * x1 * x2^2 + 1/3 * x2^3));
        
        indicesF = [elem(1) elem(2)];
        indicesF = double(indicesF);
        
        F(indicesF) += [F1 F2]';
    end
end
