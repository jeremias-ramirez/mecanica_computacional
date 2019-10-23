function [T] = esquemaTemporal(M, F, TI, et, N, model, xnode)
	dx = model.dx;	
	T = TI;	
	dt = 1;
	rhoCp = model.rhoCp;
	I = diag( ones(N, 1));
	tol = model.tol;
	e = 1;
	iterMax = model.iterMax; 
	
	switch (et)
		case 1
			fd = 1;
			if model.v ~= 0
				dt = min([dx / model.v, 0.5 * dx^2  / model.k * fd]);
			else
				dt =  0.5 * dx^2  / model.k * fd;
			end

			dt_rhoCp = dt / rhoCp;
			for n = 1:iterMax
				TN = dt_rhoCp * F + (I - dt_rhoCp * M) * TI;
				err = norm(TN - TI, 2) / norm(TN, 2);
				T = [T, TN];	
				if err < tol
					break;
				end
				TI = TN;
			end

		case 2

			rhoCp_dt = rhoCp / dt;
			k = rhoCp_dt * I + M;

			for n = 1:iterMax

				FF = F + rhoCp_dt * TI;
				TN = k \ FF;
				err = norm(TN - TI, 2) / norm(TN, 2);
				T = [T, TN];	

				if err < tol
					break;
				end
				TI = TN;
				
				
			end
	end
end

