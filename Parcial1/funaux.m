function k = getMatrixTriang(N,a,b,c)
	k = diag(b*ones(N-2,1))+diag(c*ones(N-3,1),1)+diag(a*ones(N-3,1),-1);
end

function k = getMatrixTriang2(N,a,b,c)
	k = diag(b*ones(N,1))+diag(c*ones(N-1,1),1)+diag(a*ones(N-1,1),-1);
end

function [m,b] = getCondicionBordes(m, dx, model, cb)
	b = zeros(length(m(:,1)),1);

	[m, b] = getCondicionIzq(m, b, dx, model, cb(1,:)');
	[m, b] = getCondicionDer(m, b, dx, model, cb(2,:)');
end

function k = getMatrixN(N,a,b,c)
	k = diag(b * ones(N,1)) + diag(c * ones(N-1,1),1) + diag(a * ones(N-1,1),-1);
	k(1,:) = [1, zeros(1, N-1)];
	k(N,:) = [zeros(1, N-1), 1];

end

function [M, F] = getCondicionIzq(M, F, dx, model, cbI)
	
	switch( cbI(1) )
		case 1
			F(1,1) = cbI(2);
		case 2
			F(1,1) = F(1,1) - 2 * cbI(2) / dx - model.v * cbI(2) / model.k ;
			M(1, 1) = 2 * model.k / dx^2 + model.c;
			M(1, 2) = - 2 * model.k / dx^2;

		case 3
			F(1,1) = F(1,1) + 2 * cbI(2) * cbI(3) / dx + model.v * cbI(2) * cbI(3) / model.k ;
			M(1,1) =  M(1, 1) + 2 * cbI(2) / dx + model.v * cbI(2) / model.k;
			M(1,2) = -2 * model.k / dx^2;
	end
end

			
function [M, F] = getCondicionDer(M, F, dx, model, cbD)
	N = length(F);
	
	switch( cbD(1) )
		case 1
			F(N,1) = cbD(2);
		case 2
			F(N, 1) = F(1,1) - 2 * cbD(2) / dx + model.v * cbD(2) / model.k ;
			M(N, N) = 2 * model.k / dx^2 + model.c;
			M(N, N-1) = - 2 * model.k / dx^2;
		case 3
			F(N,1) = F(1,1) + 2 * cbD(2) * cbD(3) / dx + model.v * cbD(2) * cbD(3) / model.k ;
			M(N,N) =  M(N, N) + 2 * cbD(2) / dx + model.v * cbD(2) / model.k;
			M(N,N-1) = -2 * model.k / dx^2;
	end
end



function [M, F] =  getSystem(N, dx, model, cb)
	k = model.k;
	v = model.v;
	c = model.c;
	G = model.G;
	
	% numero de Peclet
    	Pe = (v * dx)/(2*k);
   	 
    	% calculo segun el valor del numero de peclet
    	if (Pe > 1)
    	    knum = v * dx / 2;
    	else
    	    knum = 0;
    	end
	%k = k + knum ;
	
	a = -1 * (k / dx^2 + v / (2 * dx)) ;
	b = 2 * k / dx^2 + c ; 
	c = - k / dx^2 + v / (2 * dx) ;
	

	if length(G) == 1
		G = ones(N,1) * G;
	end

	F = G;

	M = getMatrixN(N, a, b, c);
	
	[M, F] = getCondicionIzq(M, F, dx, model, cb(1,:)');
	[M, F] = getCondicionDer(M, F, dx, model, cb(2,:)');
	
end

	
function [T] = showNoEstacionario(M, F, TI, et, dx, model, xnode)
	
	N = length(F);

	fd = 1 ;
	dt = 0.0;
	if model.v ~= 0
		dt = min([dx / model.v, 0.5 * dx^2  / model.k * fd]);
	else
		dt =  0.5 * dx^2  / model.k * fd;
	end
	
	if model.dt > 0 
		dt = model.dt;
	end

	rhoCp = model.rhoCp;
	
	I = diag(ones(N,1));
	
	tol = 0.00001;
	e = 1;
	iter = 0;
	iterMax = model.iterMax; 

	switch (et)
		case 1
			dt_rhoCp = dt / rhoCp;
			for n = 1:iterMax
				TN = dt_rhoCp * F + (I - dt_rhoCp * M) * TI
				err = norm(TN - TI, 2) / norm(TN, 2);
				
				if err < tol
					break;
				end
				plot(xnode, TN)
				pause(0.5)
				TI = TN;
			end

		case 2

			dt = dt;
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
				TI = TN
				
				
			end
	end
end


function k = getMatrixExpli(N,a,b,c)
	k = diag(-b*ones(N,1))+diag(c*ones(N-1,1),1)+diag(a*ones(N-1,1),-1);
	k(1,:) = [1,zeros(1,N-1)];
	k(N,:) = [zeros(1,N-1),1];
end

function k = getMatrixInter(N,a)
	k = diag(-a*ones(N,1))+diag(ones(N-1,1),1)+diag(ones(N-1,1),-1);
end

#pos true = inicio
#pos false = fin

function [m, b] = setContornoDirichlet(m,b,N,pos,T)
	if pos
	       	m(1,:) = [1,zeros(1,N-1)];
		b(1,1) = T;
	else 
		m(N,:) = [zeros(1,N-1),1];
		b(N,1) = T;
	end
end

function [m, b] = setContornoNeumann_N(m,b,N,q,k,h,orden,g)
	if orden == 2 
		#utilizando nodo ficticio eliminado
		m(N,:) = [zeros(1,N-2),-1,1];
		b(N,1) = (q*h)/k+(h^2*g)/2*k;

		#utilizando decentrada
		#m(N,:) = [zeros(1,N-3),1,-4,3];
		#b(N,1) = (q*2*h)/k;
	else 
		m(N,:) = [zeros(1,N-2),-1,1];
		b(N,1) = (q*h)/k;
	end

end

function b = getVectorInter(N,h,G,k)
	b= 1/k .* -G .* h^2 .* ones(N,1);
	b(1,1)=0;
	b(N,1)=0;
end

