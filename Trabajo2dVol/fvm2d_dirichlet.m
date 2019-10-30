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

for n = 1 : size(DIR, 1) 
	P = DIR(n, 1);
	phi_P = DIR(n, 2);
	cara = DIR(n, 3);
	R = 0;
	switch (cara)
		case 1 %Sur
			R = cells(P).ks * cells(P).as / cells(P).ds;
		case 2 %Este
			R = cells(P).ke * cells(P).ae / cells(P).de;
		case 3 %Norte
			R = cells(P).kn * cells(P).an / cells(P).dn;
		case 4 %Oeste
			R = cells(P).kw * cells(P).aw / cells(P).dw;
	end

	K(P, P) = K(P, P) + R;
	F(P) = F(P) + R * phi_P;
end



end
