function [xnode, T] = getBordesIzqVol(T, xnode, model, cbI)
	dx = model.dx;
	xnode = [model.xI xnode']';
	phi_borde = 0; 	
	switch( cbI(1) )
		case 1
			phi_borde = cbI(2);

		case 2
			phi_borde = T(1) - cbI(2) * dx / (2 * model.k);
		case 3
			den = model.k + (dx/2 * cbI(2));
			a = (cbI(2) * cbI(3) * (dx/2) ) / den;
			b = model.k / den;
			phi_borde = a + b * T(1); 

	end
	T = [phi_borde; T];
end


