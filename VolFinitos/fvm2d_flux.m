function [Q] = fvm2d_flux(PHI,cells,neighb)
% Descripción: módulo calcular el flujo de calor en todo el dominio. Se aplica la 
% Ley de Fourier y se evalúa como fluye el calor en todos los centros de celdas del
% dominio.

% Entrada:
% * PHI: vector solución. Cada elemento del vector representa un valor escalar }
% asociado a cada celda de la malla, y su posición dentro del vector depende de 
% cómo se especificó en icone.
% * neighb: matriz de vecindad.
% * cells: vector de celdas.

% Salida:
% * Q: vector de flujo de calor. Se forma de una componente en sentido x, Qx, y 
% una componente en sentido y, Qy.
% ----------------------------------------------------------------------    
    Q = [];
end
