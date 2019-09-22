function [K,F] = fdm2d_robin(K,F,xnode,neighb,ROB)
    M = size(ROB, 1);
    for n = 1 : M
        P = ROB(n, 1);
	h = ROB(n, 2);
	tInf = ROB(n, 3);
	dir = ROB(n, 4);

        % 1) para cada nodo con condicion robin sumar en la diagonal de la fila
        % P de la matriz K el valor 2*h*dist, donde dist es la distancia del %%%% mal comentado es 2 * h / dist
        % nodo P al nodo interior respecto a la normal.
	dist = 0.0;
			
	switch (dir)
		case 1
			N = neighb(P, 3);
            		dist = 1 * abs(xnode(N, 2) - xnode(P, 2));
		case 2 
			W = neighb(P, 4);
			dist = 1 * abs(xnode(W, 1) - xnode(P, 1));
		case 3
			S = neighb(P, 1);
            		dist = 1 * abs(xnode(S,2) - xnode(P,2));
		case 4
			E = neighb(P, 2);
			dist = 1 * abs(xnode(E, 1) - xnode(P, 1));

	end

	K(P, P) = K(P, P) + 2 * h / dist;

        
        % 2) para cada nodo con condicion robin sumar en el vector F el valor
        % 2*h*phi_inf/dist.

	F(P) = F(P) + 2 * h * tInf / dist;

    end
end  
