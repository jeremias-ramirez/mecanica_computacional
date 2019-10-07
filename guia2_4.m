source "volFinitos.m"

model = struct();
model.k = 1;
model.v = 0;
model.c = 0;
model.rhoCp = 1;

dx = 0.05;

cb = [[1,1,-1] ; [1,0,-1]];

xnode = [0:dx:1]';

xnodeV = [0+dx/2:dx:1-dx/2]';
model.G = 100 * xnodeV ;

model.TI = zeros(length(xnodeV),1);

et = 2;

volFinitos(xnodeV, model, cb, et);


%plot(xnodeV, T2, "*")



