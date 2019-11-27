function [localK] = fem1d_genK(nodes,young,area)
% Descripción: módulo para calcular y evaluar de forma numérica la matriz de rigidez K.
% Se utilizan funciones de forma en coordenadas naturales y se resuelve la integral
% de forma numérica utilizando cuadratura de Gauss.

% Entrada:
% * nodes: nodos (x,y) del elemento. Los elementos admisibles son de 3 o 4 nodos.
% * D: matriz constitutiva.
% * t: espesor de la placa.

% Salida:
% * localK: matriz de rigidez del elemento (local).
% ----------------------------------------------------------------------

    % METER E en algún lado
    
    E = young;
    A = area;
    l = nodes(2) - nodes(1);
    
    localK = E*A/l * [1 -1; -1 1];
end