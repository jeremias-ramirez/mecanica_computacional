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
	
	disp("ASF")
        % 1) de la estructura neighb obtener los cuatro vecinos de cada
        % nodo P.
        S = neighb(P,1);
        E = neighb(P,2);
        N = neighb(P,3);
        W = neighb(P,4);

	disp("ASF")
        % 2) obtener la distancia (a partir de la esctructura xnode) de 
        % cada nodo P hacia sus cuatro vecinos. 
        % Si es un nodo de borde (es decir que neighb es -1 en alguna de
        % sus columnas) no definir esa distncia. Luego en el stencil 
        % considerar el reemplazo del nodo ficticio tomando el stencil 
        % clasico para la derivada segunda en esa direccion.
        Rs = 0; Re = 0; Rn = 0; Rw = 0;
	dx = cells(P).dx/2;
	dy = cells(=).dy/2;
       	K(P, P) = c(P) * cells(P).v; 
	disp("ASF")
	disp("ASF")
	if(S ~= -1)
		Rs = (cells(P).as * cells(P).ks) / (cells(P).ds + dy);
		K(P, P) = K(P, P) + Rs; 
		K(P, S) = -Rs;
        end
	
	if(N ~= -1)
		Rn = (cells(P).an * cells(P).kn) / (cells(P).dn + dy);
		K(P, P) = K(P, P) + Rn; 
	    	K(P, N) = -Rn;
        end

        if(E ~= -1)
		Re = (cells(P).ae * cells(P).ke) / (cells(P).de + dx);
		K(P, P) = K(P, P) + Re; 
		K(P, E) = -Re;
        end
	

        if(W ~= -1)
		Rw = (cells(P).aw * cells(P).kw) / (cells(P).dw + dx);
		K(P, P) = K(P, P) + Rw; 
		K(P, W) = -Rw;
        end
       
        F(P) = G(P) * cells(P).v;
        
    end


end
