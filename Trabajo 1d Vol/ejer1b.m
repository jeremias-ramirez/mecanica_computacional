source "volFinitos.m"

model = struct();
model.k = 1;
model.c = 0;
model.rhoCp = 0;
model.xI = 0;
model.xF = 1;
model.G = 0 ;
model.upwind = 0;
Pe = 10;
model.v = Pe;
et = 0;
cb = [[1,1,-1] ; [1,0,-1]];

solAnalitica = @(x, Pe) (exp(Pe) - exp(Pe .* x)) / (exp(Pe) - 1);
model2 = model;

model2.upwind = 1;
Pe2 = 10;
model2.v = Pe2;

iter = 100:-1:5;
e = zeros(size(iter),1);
e2 = zeros(size(iter),1);
dxV = zeros(size(iter),1);
for i = iter 
	dx = (model.xF - model.xI)/i;
	model.dx = dx;
	model2.dx = dx;

	dxV(i-4) = dx;

	xnode = [0; [0+dx/2:dx:1-dx/2]'; 1];
	T1 = solAnalitica(xnode, Pe);
	T3 = solAnalitica(xnode, Pe2);
	[xnodeV, T2] = volFinitos(model, cb, et);
	[xnodeV, T4] = volFinitos(model2, cb, et);
	e2(i-4) = sum((T4-T3).^2) / length(xnodeV);
	e(i-4) = sum((T2-T1).^2) / length(xnodeV);
end
figure(1)
plot(log(dxV), log(e), log(dxV), log(e2))
title("Grafica log(error)")
xlabel("cantidad de celdas (N)")
ylabel("log(e)")
%saveas(1, "ejer1bpe1.jpg")

(log(e(96))-log(e(1)))/(log(dxV(96))-log(dxV(1)))

log(e2(96)-e2(1))/(log(dxV(96)-dxV(1)))

log(e2(96)) / log(dxV(96))

