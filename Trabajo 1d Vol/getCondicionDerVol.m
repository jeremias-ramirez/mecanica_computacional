function [M, F] = getCondicionDerVol(M, F, S, model, cbD)
	dx = model.dx;
	N = length(F);
	switch( cbD(1) )
		case 1
			if model.upwind == 0 
				M(N, N) = M(N, N) - model.v * S / 2 + model.k * S / dx; 		
				F(N,1) = F(N, 1) + (- model.v * S + 2 * model.k * S / dx) * cbD(2);
			else
				if model.v > 0
					M(N, N) = M(N, N) +  model.k * S / dx; 		
					F(N,1) = F(N, 1) + ( 2 * model.k * S / dx) * cbD(2);
				end
			end
		case 2
			M(N, N) = M(N, N) +   model.v * S / 2 - model.k * S / dx; 		
			F(N, 1) = F(N, 1) + (model.v * S * dx / (2 * model.k) - S ) * cbD(2);
		case 3
			den = model.k + (dx/2 * cbD(2));
			a = (cbD(2) * cbD(3) * (dx/2) ) / den;
			b = model.k / den;
			c = (cbD(3) * cbD(2)) / den; 
			d = -cbD(2) / den;

			if model.upwind == 0
								
				M(N, N) = M(N, N) -  model.v * S / 2 - model.k * S / dx + model.v * S * b - model.k * d * S; 		
				F(N,1) = F(N, 1) + model.k * c * S - model.v * S * a;
			else 
				if model.v > 0
					M(N, N) = M(N, N) - model.k * S / dx - model.v * S + model.v * S * b - model.k * d * S; 		
					F(N,1) = F(N, 1) + model.k * c * S - model.v * S * a;


				end
			end


	end
end


