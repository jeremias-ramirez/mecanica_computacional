function [K,F] = fdm2d_dirichlet(K,F,DIR)
    M = size(DIR, 1);
    for n = 1 : M
	P = DIR(n, 1);
	phi_P = DIR(n, 2);

    	% 1) Rellenar la fila DIR(n) de K con 0 y en la columna central de dicha
    	% fila colocar un uno.
  	K(P, :) = 0;
       	K(P, P) = 1;

    	% 2) Asignar en F(DIR(n)) el valor de la condicion Dirichlet 
    	F(P) = phi_P;

    end
end
