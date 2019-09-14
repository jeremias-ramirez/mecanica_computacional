function [p,b,m] = resolverDiriNeuma(N)
	xl = 1;
	xr = 2;
	k = 1;
	G = 1;
	c = 1;
	T0 = 0;
	q = 0;
	h=(xr-xl)/N;
	a = 2+c/k * h^2;
	m = getMatrixInter(N,a);
	b = getVectorInter(N,h,G,k);

	[m, b] = setContornoDirichlet(m,b,N,true,T0);
	[m1, b1] = setContornoNeumann_N(m,b,N,q,k,h,1,G);
	[m2, b2] = setContornoNeumann_N(m,b,N,q,k,h,2,G);
	p1 = m1\b1;
	p2 = m2\b2;
	x = linspace(xl,xr,N);
	
	#y = @(x) -5/2 * x .^2 + 3/2 * x + 1;
	#yd = y(x);

	plot(x, p1,x,p2)
	
end
resolverDiriNeuma(500)

