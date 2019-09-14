
function [p,b,m] = resolverDiriDiri(N)
	xl = 0;
	xr = 1;
	k = 2;
	G = 10;
	T0 = 1;
	TN = 0;
	h=(xr-xl)/N;
	
	m = getMatrixInter(N,2);
	b = getVectorInter(N,h,G,k);

	[m, b] = setContornoDirichlet(m,b,N,true,T0);
	[m, b] = setContornoDirichlet(m,b,N,false,TN);
	p = m\b;
	y = @(x) -5/2 * x .^2 + 3/2 * x + 1;
	x = linspace(xl,xr,N);
	yd = y(x);
	plot(x, p, x,yd)
	
end

resolverDiriDiri(100)
