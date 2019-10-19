function [PHI,Q] = fvm2d_explicit(K,F,cells,neighb,model,dt)
% Descripción: módulo para resolver el sistema lineal de ecuaciones utilizando 
% esquema temporal explícito.

% Entrada:
% * K: matriz del sistema (difusión + reacción)
% * F: vector de flujo térmico.
% * cells: vector de celdas.
% * neighb: matriz de vecindad.
% * model: struct con todos los datos del modelo (constantes, esquema numérico, etc.)
% * dt: paso temporal.

% Salida:
% * PHI: vector solución. Cada elemento del vector representa un valor escalar asociado
% al centroide de cada celda de la malla, y su posición dentro del vector depende de 
% cómo se especificó cada celda en icone. Se devuelve un resultado por cada iteración 
% del método (nit columnas).
% * Q: vector de flujo de calor. Se forma de una componente en sentido x, Qx, y una 
% componente en sentido y, Qy. Se devuelve un resultado por cada iteración del 
% método (2xnit columnas).
% ----------------------------------------------------------------------
I = eye(model.nnodes,model.nnodes);
    
PHI = model.PHI_n;
PHI_n = model.PHI_n;
PHI_vec = PHI;
Q_vec = zeros(model.nnodes,2);

for n = 1 : model.maxit
    PHI = A*F + (I - A*K)*PHI_n;

    err = norm(PHI-PHI_n,2)/norm(PHI,2);
    
    PHI_n = PHI;
    PHI_vec = [PHI_vec PHI];
    [Q] = fdm2d_flux(PHI,neighb,xnode,model.k);
    Q_vec = [Q_vec, Q];
    
    if err < model.tol
        disp('Método terminado por tolerancia de error.');
        return;
    end
end

disp('Método terminado por límite de iteraciones.');

PHI = [];
Q = [];
end
