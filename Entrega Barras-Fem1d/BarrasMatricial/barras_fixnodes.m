function [Kd,Fd] = barras_fixnodes(K,F,Fixnodes)
    
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

    Kd = K;
    Fd = F;

    for ii = 1 : size(Fixnodes, 1)
       
        nodo = Fixnodes(ii, 1);
        direccion = Fixnodes(ii, 2);
        desplazamiento = Fixnodes(ii, 3);
        indiceGlobal = 2*nodo + (direccion - 2);
        
        Kd(indiceGlobal, :) = 0;
        Kd(indiceGlobal, indiceGlobal) = 1;
        Fd(indiceGlobal) = desplazamiento;
    end
end
