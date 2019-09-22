function k = getMatrixTriang(N,a,b,c)
	k = diag(b*ones(N-2,1))+diag(c*ones(N-3,1),1)+diag(a*ones(N-3,1),-1);
end

function k = getMatrixN(N,a,b,c)
	k = diag(b*ones(N,1))+diag(c*ones(N-1,1),1)+diag(a*ones(N-1,1),-1);
	k(1,:) = [1, zeros(1, N-1)];
	k(N,:) = [zeros(1, N-1), 1];

end

function [m, b] = getCondicionIzq(m, b, dx, model, cbI)
	
	switch( cbI(1) )
		case 1
			b(1,1) = cbI(2);
		case 2
			b(1,1) = (2 * cbI(2) * dx / model.k) * m(2,1);
			m(1,1) = m(2,2); 
			m(1,2) = 2;
		case 3
			b(1,1) = (-2 * cbI(2) * dx / model.k) * cbI(3) * m(2,1);
			m(1,1) =  (-2 * cbI(2) * dx / model.k) * m(2,1) + m(2,2);
			m(1,2) = 2;
	end
end

			
function [m, b] = getCondicionDer(m, b, dx, model, cbD)
	N = length(b);

	switch( cbD(1) )
		case 1
			b(N,1) = cbD(2);
		case 2
			b(N,1) = (2 * cbD(2) * dx / model.k) * m(2,3);
			m(N,N) = m(2,2); 
			m(N, N-1) = 2;
		case 3
			b(N, 1) = (-2 * cbD(2) * dx / model.k) * cbD(3) * m(2,3);
			m(N, N) =  (-2 * cbD(2) * dx / model.k) * m(2,3) + m(2,2);
			m(N, N-1) = 2;
	end
end



function [m,b] = getCondicionBordes(m, dx, model, cb)
	b = zeros(length(m(:,1)),1);

	[m, b] = getCondicionIzq(m, b, dx, model, cb(1,:)');
	[m, b] = getCondicionDer(m, b, dx, model, cb(2,:)');
end
function [T] = showNoEstacionario(m, b, G, TI, et, dx, model, xnode)
	
	N = length(b);
	fd = 1 ;
	dt = 0.0;
	if model.v ~0
		dt = min([dx / model.v, 0.5 * dx^2  / model.k * fd]);
	else
		dt =  0.5 * dx^2  / model.k * fd;

	rhoCp = model.rhoCp;
	
	
	I = diag(ones(N,1));

		
	tol = 0.00001;
	e = 1;
	iter = 0;
	iterMax = 1000;

	switch (et)
	case 1
		
		dt_rhoCp = dt / rhoCp;
		k_dx2 = model.k / dx^2;
		m(2:N-1,:) = dt_rhoCp * k_dx2 * m(2:N-1,:) + I(2:N-1,:);
		G = -dt_rhoCp * k_dx2 * G;
		T =  G + m * TI + b;
		while ( e > tol || iter < iterMax)	
			TNew =  G + m * T;

			e = norm(TNew - T, 2);
			T = TNew;
			plot(xnode, T)
			pause(0.05)
		end
	case 2

		fa = 4 ;
		dt = dt * fa;
		rhoCp_dt = rhoCp / dt;
		m(2:N-1,:) = m(2:N-1,:) - rhoCp_dt * I(2:N-1,:);
		
		T =  TI;

		while ( e > tol || iter < iterMax)	
			T(1,1) = 0;
			T(N,1) = 0;

			TNew = m \ (G  - rhoCp_dt * T + b);
			TNew(1,1) = b(1,1);
			TNew(N,1) = b(N,1);

			e = norm(TNew - T, 2);
			T = TNew;
			plot(xnode, T)
			pause(0.05)
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

