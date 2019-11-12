function [F] = fvm2d_neumann(F,cells,NEU)
% Descripción: módulo para calcular y ensamblar las contribuciones de
% celdas con caras pertenecientes a fronteras de tipo Neumann.

% Entrada:
% * F: vector de flujo térmico.
% * cells: vector de celdas.
% * NEU: matriz con la información sobre la frontera de tipo Neumann. 
%     - Columna 1: índice de la celda donde se aplica la condición de borde.
%     - Columna 2: valor de flujo térmico (q) asociado al lado del elemento.
%     - Columna 3: dirección y sentido del flujo:
%         1) Flujo en dirección eje-y, sentido negativo (S – South - Sur)
%         2) Flujo en dirección eje-x, sentido positivo (E – East - Este)
%         3) Flujo en dirección eje-y, sentido positivo (N – North – Norte)
%         4) Flujo en dirección eje-x, sentido negativo (W – West – Oeste)

% Salida:
% * F: vector de flujo térmico con modificaciones luego de aplicar la 
% condición de borde.
% ----------------------------------------------------------------------


for n = 1 : size(NEU, 1) 
	P = NEU(n, 1);
	q_P = NEU(n, 2);
	dir = NEU(n, 3);

	switch (dir)
		case 1 %Sur
			F(P) = F(P) - q_P * cells(P).as;
		case 2 %Este
			F(P) = F(P) - q_P * cells(P).ae;
		case 3 %Norte
			F(P) = F(P) - q_P * cells(P).an;
		case 4 %Oeste
			F(P) = F(P) - q_P * cells(P).aw;
	end

end


end
