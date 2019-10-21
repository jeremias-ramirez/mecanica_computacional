source "volFinitos.m"

model = struct();
model.k = 1;
model.c = 1;
model.v = 20;
model.rhoCp = 0;
model.dx = 0.01;
model.xI = 0;
model.xF = 4;
model.G = 1273 ;
et = 0
cb = [[1,373,-1] ; [3,10,298]];
cb'

figure(1)
[xnodeV, T2] = volFinitos(model, cb, et);
plot(xnodeV, T2)


