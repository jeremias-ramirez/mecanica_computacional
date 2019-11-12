function [M, F] = getCondicionIzqVol(M, F, S, model, cbI)
	dx = model.dx;
	
	switch( cbI(1) )
		case 1
			if model.upwind == 0	
				M(1, 1) = M(1, 1) + model.v * S / 2 + model.k * S / dx; 		
				F(1, 1) = F(1, 1) +  (model.v * S + 2 * model.k * S / dx) * cbI(2);
			else 
				if model.v > 0
					M(1, 1) = M(1, 1) + model.k * S / dx; 
					F(1, 1) = F(1, 1) + (model.v * S + 2 * model.k * S / dx) * cbI(2);
				end
			end
		case 2
			M(1, 1) = M(1, 1) -  model.v * S / 2 - model.k * S / dx; 		
			F(1,1) = F(1, 1) + (-model.v * S * dx /  (2 * model.k) - S ) * cbI(2);
		case 3
			den = model.k + (dx/2 * cbI(2));
			a = (cbI(2) * cbI(3) * (dx/2) ) / den;
			b = model.k / den;
			c = (cbI(3) * cbI(2)) / den; 
			d = -cbI(2) / den;
			if model.upwind == 0

				M(1, 1) = M(1, 1) +  model.v * S / 2 - model.k * S / dx - model.v * S * b - model.k * d * S; 		
				F(1,1) = F(1, 1) + model.k * c * S + model.v * S * a;
			else
				if model.v > 0
					M(1, 1) = M(1, 1) - model.k * S / dx  - model.v * S * b - model.k * d * S; 		
					F(1,1) = F(1, 1) + model.k * c * S + model.v * S * a;

				end
			end

	end
end


