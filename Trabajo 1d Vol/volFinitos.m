source "getSystemVol.m"
source "getBordesIzqVol.m"
source "getBordesDerVol.m"


function [xnode, T] = volFinitos (model, cb, et)
	
	dx = model.dx;
	XI = model.xI;
	XF = model.xF;

	xnode = [XI + dx/2 : dx : XF - dx/2]';
 	
	N = length(xnode);
	
	[M, F] = getSystemVol(N, model, cb);

	if et == 0
		T = M \ F;
		[xnode, T] = getBordesIzqVol(T, xnode, model, cb(1, :)');  
		[xnode, T] = getBordesDerVol(T, xnode, model, cb(2, :)');  

		return	
	end
	
	%TI = model.TI;
	TI = zeros(N, 1);

	if cb(1,1) == 1
		TI(1,1) = cb(1,2);
	end

	if cb(2,1) == 1
		TI(N,1) = cb(2,2);
	end

	
	T = esquemaTemporal(M, F, TI, et, N, model, xnode);

end


