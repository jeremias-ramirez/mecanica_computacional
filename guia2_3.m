source "volFinitos.m"

model = struct();
model.k = 1;
model.v = 0;
model.c = 0;
model.rhoCp = 1;

dx = 0.005;

cb = [[1,0,-1] ; [1,0,-1]];

xnode = [1:dx:2]';

xnodeV = [1+dx/2:dx:2-dx/2]';
model.G = 100 ./ xnodeV .^3 ;


et = 0;

T2 = volFinitos(xnode, model, cb, et);


plot(xnodeV, T2, "*")



