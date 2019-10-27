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

	N = size(cells, 2);
    Qx = zeros(N, 1);
    Qy = zeros(N, 1);
  
    for P = 1 : N
        S = neighb(P, 1);
        E = neighb(P, 2);
        N = neighb(P, 3);
        W = neighb(P, 4);
        
        if (E ~= -1 && W ~= -1)
            dxe = cells(P).de;
            dxw = cells(P).dw;
            Qx(P) = (PHI(E) - PHI(W)) / (dxe + dxw);
        elseif (E ~= -1)
            dx = cells(P).de;
            Qx(P) = (PHI(E) - PHI(P)) / dx;
        else
            dx = cells(P).dw;
            Qx(P) = (PHI(P) - PHI(W)) / dx;
        end
        
        if (N ~= -1 && S ~= -1)
            dyn = cells(P).dn;
            dys = cells(P).ds;
            Qy(P) = (PHI(N) - PHI(S)) / (dyn + dys);
        elseif (N ~= -1)
            dy = cells(P).dn;
            Qy(P) = (PHI(N) - PHI(P)) / dy;
        else
            dy = cells(P).ds;
            Qy(P) = (PHI(P) - PHI(S)) / dy;
        end
    end
    
    Q = [Qx, Qy];

end
