function [T] = getNoEstacionarioV(M, F, TI, et, model)
	
	N = length(F);
	dx = model.dx
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
				TI = TN;
				
				
			end
	end
end

