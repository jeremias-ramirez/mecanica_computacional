source "getSystemVol.m"

function [xnode, T] = volFinitos (model, cb, et)
	
	dx = model.dx;
	XI = model.xI;
	XF = model.xF;

	xnode = [XI + dx/2 : dx : XF - dx/2]';
 	
	N = length(xnode);
	
	[M, F] = getSystemVol(N, model, cb);

	if et == 0
		T = M \ F;
		return	
	end
	
	TI = model.TI;

	if cb(1,1) == 1
		TI(1,1) = cb(1,2);
	end

	if cb(2,1) == 1
		TI(N,1) = cb(2,2);
	end

	figure(1)
	
	showNoEstacionarioV(M, F, TI, et, dx, model, xnode)

end


