function [K,F] = fem1d_fixnodes(K,F,Fixnodes)
% Descripción: módulo para calcular y ensamblar las contribuciones de nodos 
% pertenecientes a fronteras con desplazamientos fijos.

% Entrada:
% * K: matriz del sistema.
% * F: vector de fuerzas.
% * Fixnodes: matriz con la información sobre la frontera de tipo Dirchlet.
%   - Columna 1: índice de nodo
%   - Columna 2: valor del desplazamiento.

% Salida:
% * K: matriz del sistema luego de realizar las simplificaciones que surgen de 
% aplicar la condición de borde.
% * F: vector de fuerzs luego de realizar las simplificaciones que surgen de 
% aplicar la condición de borde.
% ----------------------------------------------------------------------

    for ii = 1 : size(Fixnodes, 1)
        indiceNodo = Fixnodes(ii, 1);
        desplazamiento = Fixnodes(ii, 2);
        
        K(indiceNodo, :) = 0;
        K(indiceNodo, indiceNodo) = 1;
        F(indiceNodo) = desplazamiento;
    end
end