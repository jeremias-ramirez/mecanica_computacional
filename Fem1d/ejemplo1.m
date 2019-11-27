clear all

syms E A L b P;

xnode = [
    0;
    L
];

icone = [
    1 2
];

model.nnodes = length(xnode);
model.nelem = size(icone, 1);
model.young = E;
model.area = A;

Sideload = [
    1 2 b 0
];

Pointload = [
    2 P
];

Fixnodes = [
    1 0
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