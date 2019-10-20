function [M, F] = getCondicionIzqVol(M, F, S, model, cbI)
	dx = model.dx;
	
	switch( cbI(1) )
		case 1
			M(1, 1) = M(1, 1) + model.v * S / 2 + model.k * S / dx; 		
			F(1,1) = F(1, 1) +  (model.v * S + 2 * model.k * S / dx) * cbI(2);
		case 2
			M(1, 1) = M(1, 1) -  model.v * S / 2 - model.k * S / dx; 		
			F(1,1) = F(1, 1) + (model.v * S * dx /  (2 * model.k) + S ) * cbI(2);
		case 3
			den = model.k - (dx/2 * cbI(2));
			a = (cbI(3) * cbI(2)) / den; 
			b = - cbI(2) / den;
			c = (- cbI(2) * cbI(3) * (dx/2) ) / den;
			d = model.k / den;

			M(1, 1) = M(1, 1) +  model.v * S / 2 - model.k * S / dx - model.v * S * d - model.k * b * S; 		
			F(1,1) = F(1, 1) + model.k * a * S + model.v * s * c;


	end
end


