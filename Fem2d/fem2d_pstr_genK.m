function [localK] = fem2d_pstr_genK(nodes,D,th)
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

% Si el elemento es triangular, localK tiene que ser 3x3 y B de 3x6
% Si es un cuadrángulo, localK tiene que ser 4x4 y B de 

% localK = B' * D * B

    localK = [];
end