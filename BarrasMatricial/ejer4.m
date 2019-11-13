clear all

syms E A L P;

xnode = [
    0, 0;
    0, 1;
    1, 0;
];
xnode *= L;

icone = [1 3;
        1 2;
        3 2;
];

Fixnodes = [
    1 1 0;
    1 2 0;
    2 1 0;
    3 2 0;
];

Pointload = [
    3 -P 0;
    2 0 2*P;
];

model.nnodes = size(xnode,1);
model.nbarras = size(icone,1);

propiedades(1) = struct('E', E, 'A', A);
propiedades(2) = struct('E', E, 'A', A);
propiedades(3) = struct('E', E, 'A', sqrt(2)*A);

[K,F] = barras_initialize(model.nnodes);
K = sym(K);
F = sym(F);
[K,F] = barras_gen_system(K,F,xnode,icone,propiedades,model);
[F] = barras_pointload(F,Pointload);
[Kd,Fd] = barras_fixnodes(K,F,Fixnodes);
U = Kd\Fd;
[Def,Ten] = barras_DT(xnode,icone,propiedades, model,U);
reaction = K*U;
