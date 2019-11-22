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

% Si el elemento es triangular, localK tiene que ser 6x6 y B de 3x6
% Si es un cuadrángulo, localK tiene que ser 8x8 y B de 3x8

    formanTriangulo = @(nodos) size(nodos,1) == 3;

    if formanTriangulo(nodes)
        x = nodes(:, 1);
        y = nodes(:, 2);
        
        J = [x(2)-x(1)  y(2)-y(1);
             x(3)-x(1)  y(3)-y(1)];
        DN = [-1 1 0; -1 0 1];
        V = inv(J)*DN;
        B = [V(1,1)     0     V(1,2)     0     V(1,3)     0   ;
               0      V(2,1)    0      V(2,2)    0      V(2,3);
             V(2,1)   V(1,1)  V(2,2)   V(1,2)  V(2,3)   V(1,3)];
        
        A = 0.5;    % Area del triangulo master
        localK = B' * D * B * det(J) * A * th;
    else 
        % cuatro puntos de Gauss con peso w=1
        p = [-1 -1;
              1 -1;
              1  1;
             -1  1];
        p *= sqrt(3)/3;
        localK = zeros(8);
        
        for ii = 1 : 4
            s = p(ii,1);
            t = p(ii,2);
            DNnum = [(-1+t)/4   (1-t)/4    (1+t)/4   (-1-t)/4;
                     (-1+s)/4   (-1-s)/4   (1+s)/4   (1-s)/4];
            J = DNnum*nodes;
            V = inv(J)*DNnum;
            B = [V(1,1)     0     V(1,2)     0     V(1,3)     0     V(1,4)     0;
                   0      V(2,1)    0      V(2,2)    0      V(2,3)    0      V(2,4);
                 V(2,1)   V(1,1)  V(2,2)   V(1,2)  V(2,3)   V(1,3)  V(2,4)   V(1,4)];
            
            localK += B' * D * B * det(J) * th;
        end
    end
end