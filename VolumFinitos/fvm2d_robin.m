function [K,F] = fvm2d_robin(K,F,cells,ROB)
% Descripción: módulo para calcular y ensamblar las contribuciones de 
% nodos pertenecientes a fronteras de tipo Robin.

% Entrada:
% * K: matriz del sistema (difusión + reacción)
% * F: vector de flujo térmico.
% * cells: vector de celdas.
% * ROB: matriz con la información sobre la frontera de tipo Robin.
%     - Columna 1: índice de la celda donde se aplica la condición de borde.
%     - Columna 2: valor de coeficiente de calor (h)
%     - Columna 3: valor de temperatura de referencia (phi_inf).
%     - Columna 4: dirección y sentido del flujo:
%         1) Flujo en dirección eje-y, sentido negativo (S – South – Sur)
%         2) Flujo en dirección eje-x, sentido positivo (E – East – Este)
%         3) Flujo en dirección eje-y, sentido positivo (N – North – Norte)
%         4) Flujo en dirección eje-x, sentido negativo (W – West – Oeste)

% Salida:
% * K: matriz del sistema (difusión + reacción) con modificaciones luego
% de aplicar la condición de borde.
% * F: vector de flujo térmico con modificaciones luego de aplicar la 
% condición de borde.
% ----------------------------------------------------------------------
end
