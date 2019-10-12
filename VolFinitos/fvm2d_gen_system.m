function [K,F] = fvm2d_gen_system(K,F,neighb,cells,c,G)
% Descripción: módulo para ensamblar los términos difusivo, reactivo y fuente de
% todas las celdas de la malla, generando el stencil adecuado dependiendo de si 
% alguna de sus caras pertenece o no a la frontera.

% Entrada:
% * K: matriz del sistema (difusión + reacción).
% * F: vector de flujo térmico.
% * neighb: matriz de vecindad.
% * cells: vector de celdas.
% * c: constante del término reactivo. Es un vector que permite representar c(x,y).
% * G: fuente volumétrica. Es un vector que permite representar G(x,y).

% Salida:
% * K: matriz del sistema (difusión + reacción) con modificaciones luego del ensamble.
% * F: vector de flujo térmico con modificiaciones luego del ensamble.
% ----------------------------------------------------------------------


for P = 1 : size(cells,1)
        % 1) de la estructura neighb obtener los cuatro vecinos de cada
        % nodo P.
        S = neighb(P,1);
        E = neighb(P,2);
        N = neighb(P,3);
        W = neighb(P,4);

        % 2) obtener la distancia (a partir de la esctructura xnode) de 
        % cada nodo P hacia sus cuatro vecinos. 
        % Si es un nodo de borde (es decir que neighb es -1 en alguna de
        % sus columnas) no definir esa distncia. Luego en el stencil 
        % considerar el reemplazo del nodo ficticio tomando el stencil 
        % clasico para la derivada segunda en esa direccion.
        Rs = 0; ke = 0; kn = 0; kw = 0;

        if(S ~= -1)
            ks = abs(xnode(S,2) - xnode(P,2));
        end
	
	if(N ~= -1)
            dn = abs(xnode(N,2) - xnode(P,2));
        end

        if(E ~= -1)
            de = abs(xnode(E,1) - xnode(P,1));
        end
	

        if(W ~= -1)
            dw = abs(xnode(W,1) - xnode(P,1));
        end
        
        % 3) Considerando d2T/dx2 = ax*T(i+1,j)+bx*T(i,j)+cx*T(i-1,j) y
        % d2T/dy2 = ay*T(i,j+1)+by*T(i,j)i+cy*T(i,j-1), calcular los 
        % coeficientes asociados a cada nodo vecino en base a las distancias 
        % calculadas en el paso anterior y asi definir el stencil.
        

        if (W == -1)
            ax = 2/(de*de);
            bx = -2/(de*de);
            cx = 0;
        elseif (E == -1)
            ax = 0;
            bx = -2/(dw*dw);
            cx = 2/(dw*dw);
        else
            ax = 2/(de*(de+dw));
            bx = -2/(de*dw);
            cx = 2/(dw*(de+dw));
        end


        if (S == -1)
            ay = 2/(dn*dn);
            by = -2/(dn*dn);
            cy = 0;
        elseif (N == -1)
            ay = 0;
            by = -2/(ds*ds);
            cy = 2/(ds*ds);
        else
            ay = 2/(dn*(dn+ds));
            by = -2/(ds*dn);
            cy = 2/(ds*(ds+dn));
        end

	% la funcion anterior deberia estar separada, hacer una que solo genere los coeficientes

        % (COMPLETAR LOS COEFICIENTES ay, by y cy)
        
        % 4) Asignar en la matriz K los correspondientes valores
        % calculados previamente. Asignar en el vector F el valor de la fuente
        % en el nodo P.

	
	%if ((0.5 - xnode(P,1))< 1e-05)
	%	k(P) = 0.25;
	%end
	%
	%if (xnode(P,1) > 1e-05)
	%	k(P) = 0.25;
	%end

        K(P,P) = -k(P) * (bx + by) + c(P);
	if ( W ~= -1)
		K(P, W) = -k(P) * cx;
	end
	if ( E ~= -1)
		K(P, E) = -k(P) * ax;
	end
	
	if ( N ~= -1)
		K(P, N) = -k(P) * ay;
	end
	
	if ( S ~= -1)
		K(P, S) = -k(P) * cy;
	end
	
	%if ((xnode(P,1) - 0.5) >= 1e-05)
	%	G(P) = -500;
	%end
	%
	%if ((1.0 - xnode(P,1)) < 1e-05)
	%	G(P) = -500;
	%end


        F(P) = G(P);
        
    end


end
