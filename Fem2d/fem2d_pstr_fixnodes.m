function [K,F] = fem2d_pstr_fixnodes(K,F,Fixnodes)
% Descripción: módulo para calcular y ensamblar las contribuciones de nodos 
% pertenecientes a fronteras con desplazamientos fijos.

% Entrada:
% * K: matriz del sistema.
% * F: vector de fuerzas.
% * Fixnodes: matriz con la información sobre la frontera de tipo Dirchlet.
%   - Columna 1: índice de nodo
%   - Columna 2: dirección de desplazamiento:
%       1- sentido eje-x.
%       2- sentido eje-y.
%   - Columna 3: valor del desplazamiento.

% Salida:
% * K: matriz del sistema luego de realizar las simplificaciones que surgen de 
% aplicar la condición de borde.
% * F: vector de fuerzs luego de realizar las simplificaciones que surgen de 
% aplicar la condición de borde.
% ----------------------------------------------------------------------

    for i = 1 : size(Fixnodes,1)
        iNodoL = Fixnodes(i,1);
        iNodoG = 2*iNodoL - (2 - Fixnodes(i,2));
        desp = Fixnodes(i,3);
        
        K(iNodoG, :) = 0;
        K(iNodoG, iNodoG) = 1;
        F(iNodoG) = desp;
    end
end