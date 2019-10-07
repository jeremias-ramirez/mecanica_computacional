source "volFinitos.m"

model = struct();
model.k = 1;
model.v = 0;
model.c = 1;
model.rhoCp = 1;

cb = [[1,0,-1] ; [2,0,-1]];

xnode = [1:0.1:2]';

model.G = 1 ;
et = 0;

T2 = volFinitos(xnode, model, cb, et);

xnodeV = [1+0.1/2:0.1:2-0.1/2]';

plot(xnodeV, T2, "*")



