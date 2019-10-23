function [xnode, T] = getBordesDerVol(T, xnode, model, cbD)
	dx = model.dx;
	xnode = [xnode' model.xF]';
	phi_borde = 0;
	switch( cbD(1) )
		case 1
			phi_borde = cbD(2);

		case 2
			phi_borde = T(end) - cbD(2) * dx / (2 * model.k);
		case 3
			den = model.k + (dx/2 * cbD(2));
			a = (cbD(2) * cbD(3) * (dx/2) ) / den;
			b = model.k / den;
			phi_borde = a + b * T(end);

	end
	T = [T; phi_borde];
end


