model = struct();
model.k = 1;
model.c = 0;
model.rhoCp = 0;
model.xI = 0;
model.xF = 1;
model.G = 0 ;
model.upwind = 1;
Pe = 10;
model.v = Pe;
et = 0;
cb = [[1,1,-1] ; [1,0,-1]];

solAnalitica = @(x, Pe) (exp(Pe) - exp(Pe .* x)) / (exp(Pe) - 1);

iter = 5:100;
e = zeros(size(iter),1);

for i = iter 
	dx = (model.xF - model.xI)/i;
	model.dx = dx;

	xnode = [0; [0+dx/2:dx:1-dx/2]'; 1];
	T1 = solAnalitica(xnode, Pe);
	[xnodeV, T2] = volFinitos(model, cb, et);
	e(i-4) = norm(T2-T1,2);
end
figure(2)
plot(iter, log(e))
title("Grafica log(error) - Esquema Upwind - Pe = 10")
xlabel("cantidad de celdas (N)")
ylabel("log(e)")
saveas(1, "ejer1bpe1.jpg")

(log(e(96))-log(e(1)))/95