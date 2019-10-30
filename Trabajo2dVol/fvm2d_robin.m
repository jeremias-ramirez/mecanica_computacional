function [K,F] = fvm2d_robin(K,F,cells,ROB)
% Descripción: módulo para calcular y ensamblar las contribuciones de 
% nodos pertenecientes a fronteras de tipo Robin.

% Entrada:
% * K: matriz del sistema (difusión + reacción)
% * F: vector de flujo térmico.
% * cells: vector de celdas.
% * ROB: matriz con la información sobre la frontera de tipo Robin.
%     - Columna 1: índice de la celda donde se aplica la condición de borde.
%     - Columna 2: valor de coeficiente de calor (h)
%     - Columna 3: valor de temperatura de referencia (phi_inf).
%     - Columna 4: dirección y sentido del flujo:
%         1) Flujo en dirección eje-y, sentido negativo (S – South – Sur)
%         2) Flujo en dirección eje-x, sentido positivo (E – East – Este)
%         3) Flujo en dirección eje-y, sentido positivo (N – North – Norte)
%         4) Flujo en dirección eje-x, sentido negativo (W – West – Oeste)

% Salida:
% * K: matriz del sistema (difusión + reacción) con modificaciones luego
% de aplicar la condición de borde.
% * F: vector de flujo térmico con modificaciones luego de aplicar la 
% condición de borde.
% ----------------------------------------------------------------------
for n = 1 : size(ROB, 1) 
	
	P = ROB(n, 1);
	h_P = ROB(n, 2);
	phi_inf = ROB(n, 3);
	dir = ROB(n, 4);
	k = 0;
	f = 0;
	switch (dir)
		case 1 %Sur
			k = -cells(P).ks * cells(P).as * -(h_P / (cells(P).ks + h_P * cells(P).ds));
			f = cells(P).ks * cells(P).as * (h_P * phi_inf / (cells(P).ks + h_P * cells(P).ds));
		case 2 %Este
			k = -cells(P).ke * cells(P).ae * -h_P / (cells(P).ke + h_P * cells(P).de);
			f = cells(P).ke * cells(P).ae * h_P * phi_inf / (cells(P).ke + h_P * cells(P).de);
		case 3 %Norte
			k = -cells(P).kn * cells(P).an * -h_P / (cells(P).kn + h_P * cells(P).dn);
			f = cells(P).kn * cells(P).an * h_P * phi_inf / (cells(P).kn + h_P * cells(P).dn);
		case 4 %Oeste
			k = -cells(P).kw * cells(P).aw * -h_P / (cells(P).kw + h_P * cells(P).dw);
			f = cells(P).kw * cells(P).aw * h_P * phi_inf / (cells(P).kw + h_P * cells(P).dw);
	end
	K(P, P) = K(P, P) + k;
       	F(P) = F(P) + f;	

end
