function [D] = fem2d_pstr_const_mat(model)
% Descripción: módulo para calcular la matriz constitutiva del sistema. La misma 
% depende de los valores del Módulo de Young y el Coeficiente de Poisson. La forma
% de esta matriz es siempre de 3x3, pero la selección de coeficientes determina el
% tratamiento del sistema como Deformación o Tensión Plana

% Entrada:
% * model: struct con todos los datos del modelo (constantes, esquema numérico, etc.)

% Salida:
% * D: matriz constitutiva del sistema.
% ----------------------------------------------------------------------

    D = [];
    
end