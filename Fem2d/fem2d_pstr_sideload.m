function [F] = fem2d_pstr_sideload(F,Sideload,xnode,th)
% Descripción: módulo para calcular y ensamblar las contribuciones de pares de nodos
% (lados) donde se aplican cargas distribuidas.

% Entrada:
% * F: vector de fuerzas.
% * sideload: matriz con la información sobre fronteras con cargas distribuidas. 
%   - Columnas 1-2: dos nodos contiguos formando un lado de un elemento.
%   - Columna 3: valor de fuerza en sentido eje-x.
%   - Columna 4: valor de fuerza en sentido eje-y.
% * xnode: matriz de nodos con pares (x,y) representando las coordenadas de cada nodo
%   de la malla.
% * th: espesor de la placa.

% Salida:
% * F: vector de fuerzas. Presenta modificaciones luego de aplicar la condición de borde.
% ----------------------------------------------------------------------

% Acá no importa si el elemento es triangular o cuadrangular. Seguro?

    for i = 1 : size(Sideload, 1)
        iNodo1 = Sideload(i,1);
        iNodo2 = Sideload(i,2);
        fx = Sideload(i,3);
        fy = Sideload(i,4);
        
        long = norm(xnode(iNodo1,:) - xnode(iNodo2,:));
        
        F([(2*iNodo1-1) 2*iNodo1 (2*iNodo2-1) 2*iNodo2]) = F([(2*iNodo1-1) 2*iNodo1 (2*iNodo2-1) 2*iNodo2]) + ...
        0.5 * long * [fx fy fx fy]' * th;
    end
end
