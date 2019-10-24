source "volFinitos.m"

model = struct();
model.k = 1;
model.c = 0;
model.rhoCp = 0;
model.xI = 0;
model.xF = 1;
model.G = 0 ;
model.upwind = 0;

model.dx = (model.xF - model.xI)/5;
et = 0
cb = [[1,1,-1] ; [1,0,-1]];
cb'

solAnalitica = @(x, Pe) (exp(Pe) - exp(Pe .* x)) / (exp(Pe) - 1) 
xnode = [0:0.01:1]';


Pe = 0.1;
model.v = Pe 
T1 = solAnalitica(xnode, Pe);
[xnodeV, T2] = volFinitos(model, cb, et);

grid on

figure(1)
plot(xnode,T1, xnodeV, T2, "*")
title("Figura con Pe = 0.1")
xlabel("x")
ylabel("Temperatura")
legend("Sol. analitica", "Sol. aproximada")
%saveas(1, "ejer1pe01.jpg")

Pe = 1;
model.v = Pe 
T1 = solAnalitica(xnode, Pe);
[xnodeV, T2] = volFinitos(model, cb, et);

figure(2)
plot(xnode,T1, xnodeV, T2, "*")
title("Figura con Pe = 1")
xlabel("x")
ylabel("Temperatura")
legend("Sol. analitica", "Sol. aproximada")
%saveas(2, "ejer1pe1.jpg")

Pe = 10;
model.upwind = 1;
model.v = Pe 
T1 = solAnalitica(xnode, Pe);
[xnodeV, T2] = volFinitos(model, cb, et);

figure(3)
plot(xnode,T1, xnodeV, T2, "*")
title("Figura con Pe = 10, Upwind")
xlabel("x")
ylabel("Temperatura")
legend("Sol. analitica", "Sol. aproximada")
%saveas(3, "ejer1pe10.jpg")
