clear all

xnode = [
    0, 0;
    0.7, 0;
    1.4, 0;
    2.1, 0;
];

icone = [1 2;
        2 3;
        3 4;
];

Fixnodes = [
    1 1 0;
    1 2 0;
    4 1 0;
    4 2 0;
];

Pointload = [2 13500 0];

model.nnodes = size(xnode,1);
model.nbarras = size(icone,1);

propiedades(1) = struct('E', 20e9, 'A', 6e-4);
propiedades(2) = struct('E', 20e9, 'A', 6e-4);
propiedades(3) = struct('E', 10e9, 'A', 12e-4);

[K,F] = barras_initialize(model.nnodes);
[K,F] = barras_gen_system(K,F,xnode,icone,propiedades,model);
[F] = barras_pointload(F,Pointload);
[Kd,Fd] = barras_fixnodes(K,F,Fixnodes);
U = Kd\Fd;
[Def,Ten] = barras_DT(xnode,icone,propiedades, model,U);
reaction = K*U;
