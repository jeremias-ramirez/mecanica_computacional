function [M, F] = getCondicionDerVol(M, F, S, model, cbD)
	dx = model.dx;
	N = length(F);
	switch( cbD(1) )
		case 1
			M(N, N) = M(N, N) - model.v * S / 2 + model.k * S / dx; 		
			F(N,1) = F(N, 1) + (- model.v * S + 2 * model.k * S / dx) * cbD(2);
		case 2
			M(N, N) = M(N, N) +   model.v * S / 2 - model.k * S / dx; 		
			F(N, 1) = F(N, 1) - (model.v * S * dx / (2 * model.k) + S ) * cbD(2);
		case 3
			den = model.k - (dx/2 * cbD(2));
			a = (cbD(3) * cbD(2)) / den; 
			b = cbD(2) / den;
			c = (cbD(2) * cbD(3) * (dx/2) ) / den;
			d = model.k / den;

			M(N, N) = M(N, N) -  model.v * S / 2 - model.k * S / dx + model.v * S * d + model.k * b * S; 		
			F(N,1) = F(N, 1) + model.k * a * S + model.v * S * c;


	end
end


