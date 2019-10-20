source "volFinitos.m"

model = struct();
model.k = 1;
model.c = 0;
model.rhoCp = 0;
model.dx = 0.1;
model.xI = 0;
model.xF = 1;
model.G = 0 ;
et = 0
cb = [[1,1,-1] ; [1,0,-1]];
cb'

solAnalitica = @(x, Pe) (exp(Pe) - exp(Pe .* x)) / (exp(Pe) - 1) 
xnode = [0:0.1:1]';

figure(1)
Pe = 0.1;
model.v = Pe 
T1 = solAnalitica(xnode, Pe);
[xnodeV, T2] = volFinitos(model, cb, et);
plot(xnode,T1, xnodeV, T2, "*")

figure(2)
Pe = 1;
model.v = Pe 
T1 = solAnalitica(xnode, Pe);
[xnodeV, T2] = volFinitos(model, cb, et);
plot(xnode,T1, xnodeV, T2, "*")

figure(3)
Pe = 10;
model.v = Pe 
T1 = solAnalitica(xnode, Pe);
[xnodeV, T2] = volFinitos(model, cb, et);
plot(xnode,T1, xnodeV, T2, "*")

