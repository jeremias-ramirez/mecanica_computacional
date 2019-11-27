clear all

syms E A Kx L Q;

xnode = [
    0;
    L/2;
    L;
    3*L/2;
    2*L;
    5*L/2;
    3*L
];

icone = [
    1 2;
    2 3;
    3 4;
    4 5;
    5 6;
    6 7;
];

model.nnodes = length(xnode);
model.nelem = size(icone, 1);
model.young = E;
model.area = A;

alpha = -Q / (2*L);
beta = Q / (2 * L^2);

Sideload = [
    3 4 alpha beta;
    4 5 alpha beta;
    5 6 alpha beta;
    6 7 alpha beta;
];

Pointload = [
];

Fixnodes = [
    1 0;
    7 0
];

[K,F] = fem1d_initialize(model.nnodes);
K = sym(K);
F = sym(F);
K = fem1d_gen_system(K,xnode,icone,model);
F = fem1d_sideload(F,Sideload,xnode);
F = fem1d_pointload(F,Pointload);
[Kd,Fd] = fem1d_fixnodes(K,F,Fixnodes);
U = Kd\Fd;
reaction = K * U - F;