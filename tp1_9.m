source "difFinitas.m"

model = struct();
model.k = 1;
model.v = 0;
model.c = 0;
model.rhoCp = 1;

cb = [[1,1,-1];[1,0,-1]];

xnode = [0:0.1:1]';

model.G = 100 .* xnode;
et = 2;

T = difFinitas(xnode, model, cb, et);
plot(xnode,T)
