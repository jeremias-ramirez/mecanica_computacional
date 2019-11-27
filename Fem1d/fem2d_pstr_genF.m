function [localF] = fem2d_pstr_genF(nodes,Fg,th)
% Descripción: módulo para calcular el vector de fuerzas F para cada elemento,
% producto de la presencia de una fuerzas por unidad de volumen (generalmente 
% la gravedad) en dicho elemento. La integral se resuelve mediante cuadratura
% de punto medio, y se requiere evaluar el área del elemento.

% Entrada:
% * nodes: nodos (x,y) del elemento. Los elementos admisibles son de 3 o 4 nodos.
% * Fg: magnitud de la fuerza gravitatoria.
% * th: espesor de la placa.

% Salida:
% * localF: vector local de fuerzas.
% ----------------------------------------------------------------------

% Para elementos triangulares, localF tiene que ser de 6x1
% Para cuadrángulos, tiene que ser 8x1

    formanTriangulo = @(nodos) size(nodos,1) == 3;

    if formanTriangulo(nodes)
        nodo1 = nodes(1, :);
        nodo2 = nodes(2, :);
        nodo3 = nodes(3, :);
        
        vector1 = [nodo2-nodo1 0];
        vector2 = [nodo3-nodo2 0];
        A = norm(cross(vector1, vector2)) / 2;
        
        localF = -Fg * A * th / 3 * [0 1 0 1 0 1]';
    else
        nodo1 = nodes(1, :);
        nodo2 = nodes(2, :);
        nodo3 = nodes(3, :);
        nodo4 = nodes(4, :);
        
        triang1vector1 = [nodo2-nodo1 0];
        triang1vector2 = [nodo3-nodo2 0];
        A1 = norm(cross(triang1vector1, triang1vector2)) / 2;
        triang2vector1 = [nodo3-nodo1 0];
        triang2vector2 = [nodo4-nodo3 0];
        A2 = norm(cross(triang2vector1, triang2vector2)) / 2;
        A = A1 + A2;
    
        localF = -Fg * A * th / 4 * [0 1 0 1 0 1 0 1]';
    end
end