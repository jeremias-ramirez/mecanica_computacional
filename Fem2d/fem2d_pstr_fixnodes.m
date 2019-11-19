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

    for ii = 1 : size(Fixnodes, 1)
        indiceNodo = Fixnodes(ii, 1);
        direccion = Fixnodes(ii, 2);
        desplazamiento = Fixnodes(ii, 3);
        indiceGlobal = 2*indiceNodo + (direccion - 2);
        
        K(indiceGlobal, :) = 0;
        K(indiceGlobal, indiceGlobal) = 1;
        F(indiceGlobal) = desplazamiento;
    end
end