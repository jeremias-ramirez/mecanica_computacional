source "volFinitos.m"

model = struct();
model.k = 1;
model.c = 0;
model.rhoCp = 0;
model.xI = 0;
model.xF = 1;
model.G = 0 ;
model.upwind = 0;

et = 0
cb = [[1,1,-1] ; [1,0,-1]];
cb'
dx = (model.xF - model.xI)/5;
model.dx = dx;
xnode = [0; [0+dx/2:dx:1-dx/2]'; 1];

solAnalitica = @(x, Pe) (exp(Pe) - exp(Pe .* x)) / (exp(Pe) - 1) 


Pe = 0.1;
model.v = Pe ;
T1 = solAnalitica(xnode, Pe);
[xnodeV, T2] = volFinitos(model, cb, et);
disp "Pecle 0.1"
T2

grid on

figure(1)
plot(xnode,T1, xnodeV, T2, "*")
title("Figura con Pe = 0.1")
xlabel("x")
ylabel("Temperatura")
legend("Sol. analitica", "Sol. aproximada")
saveas(1, "ejer1pe01.jpg")

Pe = 1;
model.v = Pe;
T1 = solAnalitica(xnode, Pe);
[xnodeV, T2] = volFinitos(model, cb, et);
disp "Pecle 1"
T2

figure(2)
plot(xnode,T1, xnodeV, T2, "*")
title("Figura con Pe = 1")
xlabel("x")
ylabel("Temperatura")
legend("Sol. analitica", "Sol. aproximada")
saveas(2, "ejer1pe1.jpg")

Pe = 10;
model.upwind = 1;
model.v = Pe; 
T1 = solAnalitica(xnode, Pe);
[xnodeV, T2] = volFinitos(model, cb, et);
disp "Pecle 10"
T2

figure(3)
plot(xnode,T1, xnodeV, T2, "*")
title("Figura con Pe = 10, Upwind")
xlabel("x")
ylabel("Temperatura")
legend("Sol. analitica", "Sol. aproximada")
saveas(3, "ejer1pe10.jpg")



