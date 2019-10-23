source "getMatrixTriang.m"
source "getCondicionIzqVol.m"
source "getCondicionDerVol.m"

function [M, F] =  getSystemVol(N, model, cb)
	k = model.k;
	v = model.v;
	c = model.c;
	G = model.G;
	dx = model.dx;
	dy = 1;
	t = 1;
	A = dy * t;
	V = A * model.dx;
	if model.upwind == 0
		a = -A * v / 2 - k * A / dx; 
		b = c * V +  2 * k * A / dx ;
		c = v * A / 2 - k * A / dx;
	else
		if model.v > 0
			a = -A * v - k * A / dx; 
			b = c * V +  2 * k * A / dx + model.v * A;
			c = - k * A / dx;
		end
		%falta para v < 0

	end

	if length(G) == 1
		G = ones(N,1) * G;
	end

	F = V * G;

	M = getMatrixTriang(N, a, b, c);
	
	[M, F] = getCondicionIzqVol(M, F, A, model, cb(1,:)');
	[M, F] = getCondicionDerVol(M, F, A, model, cb(2,:)');

end

