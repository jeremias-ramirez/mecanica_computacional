function k = getMatrixTriang(N,a,b,c)
	k = diag(b*ones(N-2,1))+diag(c*ones(N-3,1),1)+diag(a*ones(N-3,1),-1);
end


function [p] = resolverExpli(N)
	xl = 0
	xr = 1	
	T0 = 1
	TN = 0
	k = 1

	dx = (xr-xl)/N
	fd = 1 #???
	dt = 0.5 * dx^2 * fd
	a = k * dt / dx^2
	b = -2 * k * dt / dx^2 + 1
	c = a;

	T = zeros(N,1);
	T(1,1) = T0;
	T(N,1) = TN;

	G = dt * 100 * (linspace(xl,xr,N))';
	
	m = getMatrixTriang(N,a,b,c);
	
	tol = 0.0001;
	TNew = zeros(N, 1);
	contor = zeros(N-2, 1);
	contor(1, 1) = T0;
	contor(end, 1) = TN;

	plot(linspace(xl,xr,N),T)
	while(e > tol)
		TNew(2:N-1,1) = m * T(2:N-1,1) + G(2:N-1,1) + a*contor;
		TNew(1,1) = T0;
		TNew(N,1) = TN;

		e = norm(TNew-T,2);
		T = TNew;
		plot(linspace(xl,xr,N),T)
	end

end
figure(1)
resolverExpli(50);

function [p] = resolverImpli(N)
	xl = 0
	xr = 1	
	T0 = 1
	TN = 0
	k = 1
	
	dx = (xr-xl)/N
	fd = 1 #???
	fa = 2 #??
	dt = 0.5 * dx^2 * fd * fa
	a = -k  / dx^2
	b = 2 * k  / dx^2 + 1/dt
	c = a;

	T = zeros(N,1);
	T(1,1) = T0;
	T(N,1) = TN;

	G =  100 .* (linspace(xl,xr,N))';
	G(1,1) = 0;
	G(N,1) = 0;
	
	m = getMatrixTriang(N,a,b,c);
	
	tol = 0.0001;
	TNew = zeros(N, 1);
	contor = zeros(N-2, 1);
	contor(1, 1) = -a*T0;
	contor(end, 1) = -c*TN;

	plot(linspace(xl,xr,N),T)
	while(e > tol)
		TNew(2:N-1,1) = m \ (T(2:N-1,1)/dt + G(2:N-1,1) + contor);
		TNew(1,1) = T0;
		TNew(N,1) = TN;

		e = norm(TNew-T,2);
		T = TNew;
		plot(linspace(xl,xr,N),T)
	end

end
figure(2)
resolverImpli(50)

