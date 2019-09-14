function [p,b,m] = resolverDiriDiri(N)
	xl = 1;
	xr = 2;
	k = 1;
	G = 100*(linspace(xl,xr,N))'.^3
	T0 = 0;
	TN = 0;
	h=(xr-xl)/N;
	m = getMatrixInter(N,2);
	b = getVectorInter(N,h,G,k);

	[m, b] = setContornoDirichlet(m,b,N,true,T0);
	[m, b] = setContornoDirichlet(m,b,N,false,TN);
	[m2, b2] = setContornoNeumann_N(m,b,N,1,k,h,2,G(N,1));

	p = m\b;
	y = @(x) -5/2 * x .^2 + 3/2 * x + 1;
	x = linspace(xl,xr,N);
	yd = y(x);
	plot(x, p)
	
end

resolverDiriDiri(100)

