source "difFinitas.m"

model = struct();
model.k = 1;
model.v = 3;
model.c = 0;
model.rhoCp = 1;

cb = [[1,100,-1];[2,5,-1]];

xnode = [0:1/3:1]'
pause()

model.G = 100;
et = 2;

difFinitas(xnode, model, cb, et);
