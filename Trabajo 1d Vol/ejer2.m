source "volFinitos.m"

model = struct();
model.k = 1;
model.c = 1;
model.v = 20;
model.rhoCp = 1;
model.xI = 0;
model.xF = 4;
model.G = 1273 ;
model.tol = 1e-10;
model.iterMax = 100;
model.upwind = 1;
et = 0
cb = [[1,373,-1] ; [3,10,298]];
cb'

%model.dx = (model.xF - model.xI) / 10;
%[xnodeV1, T1] = volFinitos(model, cb, et);
%T1
%model.dx = (model.xF - model.xI) / 20;
%[xnodeV2, T2] = volFinitos(model, cb, et);
%T2
%model.dx = (model.xF - model.xI) / 50;
%[xnodeV3, T3] = volFinitos(model, cb, et);
%T3
%figure(1)
%plot(xnodeV1, T1, xnodeV2, T2, xnodeV3, T3 )
%title("Soluciones")
%xlabel("x")
%ylabel("Temperatura")
%h = legend("Sol. aproximada 10 celdas", "Sol. aproximada 20 celdas", "Sol. aproximada 50 celdas");
%legend(h, "location", "northwest");
%%saveas(1, "ejer2Funciones.jpg")

%Prueba temporal
et = 2
model.dx = (model.xF - model.xI) / 50;
[xnodeV, T] = volFinitos(model, cb, et);
hold on
for i = 1 : size(T,2)
	plot(xnodeV, T(:, i))
	pause(0.1)
end

