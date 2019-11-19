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

    esTensionPlana = @(model) model.pstrs == 1;

    E = model.young;
    v = model.poiss;

    if esTensionPlana(model)
        D = E / (1- v^2) * [1 v 0 ; v 1 0 ; 0 0 (1-v)/2];
    else
        aux1 = E*(1-v) / ((1 +v ) * (1 - 2 * v));
        aux2 = v / ( 1- v);
        aux3 = (1 - 2 * v)/ (2 * (1-v));
        D = aux1 * [1 aux2 0 ; aux2 1 0 ; 0 0 aux3];
    end
end