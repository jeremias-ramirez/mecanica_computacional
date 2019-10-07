source "difFinitas.m"
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

T1 = difFinitas(xnode, model, cb, et);
T2 = volFinitos(xnode, model, cb, et);

xnodeV = [1+0.1/2:0.1:2-0.1/2]';

plot(xnode,T1, xnodeV, T2, "*")



