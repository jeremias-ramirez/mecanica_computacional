source "difFinitas.m"

model = struct();
model.k = 1;
model.v = 3;
model.c = 0;
model.rhoCp = 1;
model.dt = 1
model.iterMax = 1

cb = [[1,100,-1];[2,5,-1]];

xnode = [0:1/3:1]';

model.TI = 100 - 50 * xnode;

model.G = 100;
et = 2;

difFinitas(xnode, model, cb, et);
