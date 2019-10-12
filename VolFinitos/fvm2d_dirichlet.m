function [K,F] = fvm2d_dirichlet(K,F,cells,DIR)
% Descripción: módulo para calcular y ensamblar las contribuciones de celdas
% pertenecientes a fronteras de tipo Dirichlet.

% Entrada:
% * K: matriz del sistema (difusión + reacción)
% * F: vector de flujo térmico.
% * cells: vector de celdas.
% * DIR: matriz con la información sobre la frontera de tipo Dirchlet.
%   - Columna 1: índice de la celda donde se aplica la condición de borde.
%   - Columna 2: valor en la cara de la celda (escalar)
%   - Columna 3: cara a la que se aplica la condición de borde:
%       1) S – South – Sur
%       2) E – East – Este
%       3) N – North – Norte
%       4) W – West – Oeste

% Salida:
% * K: matriz del sistema (difusión + reacción) con modificaciones luego de 
%   aplicar la condición de borde.
% * F: vector de flujo térmico con modificaciones luego de aplicar la condición
%   de borde.
% ----------------------------------------------------------------------
end
