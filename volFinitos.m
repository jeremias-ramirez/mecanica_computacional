
source "funaux.m"

function [T] = showNoEstacionarioV(M, F, TI, et, dx, model, xnode)
	
	N = length(F);

	fd = 1 ;
	dt = 0.0;
	if model.v ~= 0
		dt = min([dx / model.v, 0.5 * dx^2  / model.k * fd]);
	else
		dt =  0.5 * dx^2  / model.k * fd;
	end
	
	rhoCp = model.rhoCp;
	
	I = diag(ones(N,1));
	
	tol = 1.00e-8;
	e = 1;
	iter = 0;
	iterMax = 1000;
	
	switch (et)
		case 1
			dt_rhoCp = dt / rhoCp;
		
			for n = 1:iterMax
				TN = dt_rhoCp * F + (I - dt_rhoCp * M) * TI;
				err = norm(TN - TI, 2) / norm(TN, 2);
				
				if err < tol
					break;
				end
				plot(xnode, TN)
				pause(0.1)
				TI = TN;
			end

		case 2

			fa = 2 ;
			dt = dt * fa;
			rhoCp_dt = rhoCp / dt;
			k = rhoCp_dt * I + M;

			for n = 1:iterMax

				FF = F + rhoCp_dt * TI;
				TN = k \ FF;
				err = norm(TN - TI, 2) / norm(TN, 2);
				
				if err < tol
					break;
				end
				plot(xnode, TN)
				pause(0.05)
				TI = TN;
				
				
			end
	end
end


function [M, F] = getCondicionIzqVol(M, F, dx, S, model, cbI)
	
	switch( cbI(1) )
		case 1
			M(1, 1) = M(1, 1) + model.v * S / 2 + model.k * S / dx; 		
			F(1,1) = F(1, 1) - (- model.v * S - 2 * model.k * S / dx) * cbI(2);
		case 2
			M(1, 1) = M(1, 1) -  model.v * S / 2 - model.k * S / dx; 		
			F(1,1) = F(1, 1) + (model.v * S * dx / model.k - S ) * cbI(2);
		case 3
			den = model.k - (dx/2 * cbI(2)
			a = (cbI(3) * cbI(2)) / den; 
			b = - cbI(2) / den;
			c = (- cbI(2) * cbI(3) * (dx/2) ) / den;
			d = model.k / den;

			M(1, 1) = M(1, 1) +  model.v * S / 2 - model.k * S / dx - model.v * S * d - model.k * b * S; 		
			F(1,1) = F(1, 1) + model.k * a * S + model.v * s * c;


	end
end

function [M, F] = getCondicionDerVol(M, F, dx, S, model, cbD)
	N = length(F);
	switch( cbD(1) )
		case 1
			M(N, N) = M(N, N) - model.v * S / 2 + model.k * S / dx; 		
			F(N,1) = F(N, 1) - (model.v * S - 2 * model.k * S / dx) * cbD(2);
		case 2
			M(N, N) = M(N, N) +   model.v * S / 2 - model.k * S / dx; 		
			F(N, 1) = F(N, 1) - (model.v * S * dx / model.k + S ) * cbD(2);
		case 3
			den = model.k - (dx/2 * cbD(2)
			a = (cbD(3) * cbD(2)) / den; 
			b = - cbD(2) / den;
			c = (- cbD(2) * cbD(3) * (dx/2) ) / den;
			d = model.k / den;

			M(N, N) = M(N, N) -  model.v * S / 2 - model.k * S / dx + model.v * S * d - model.k * b * S; 		
			F(N,1) = F(N, 1) + model.k * a * S - model.v * s * c;


	end
end

function [M, F] =  getSystemVol(N, dx, model, cb, xnode)
	k = model.k;
	v = model.v;
	c = model.c;
	G = model.G;
	dy = 1;
	t = 1;
	S = dy * t;
	V = S * dx;
	

	a = -S * v / 2 - k * S / dx; 
	b = c * V +  2 * k * S / dx ;
	c = v * S / 2 - k * S / dx;

	if length(G) == 1
		G = ones(N,1) * G;
	end

	F = V * G;

	M = getMatrixTriang2(N, a, b, c);
	
	[M, F] = getCondicionIzqVol(M, F, dx, S, model, cb(1,:)');
	[M, F] = getCondicionDerVol(M, F, dx, S, model, cb(2,:)');

end

function [T] = volFinitos (xnode, model, cb, et)
	
	N = length(xnode);
	dx = xnode(2,1) - xnode(1,1);
	
	[M, F] = getSystemVol(N, dx, model, cb, xnode);
	if et == 0
		T = M \ F;
		return	
	end
	
	TI = model.TI;

	if cb(1,1) == 1
		TI(1,1) = cb(1,2);
	end

	if cb(2,1) == 1
		TI(N,1) = cb(2,2);
	end

	figure(1)
	
	showNoEstacionarioV(M, F, TI, et, dx, model, xnode)

end


