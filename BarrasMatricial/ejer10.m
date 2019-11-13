clear all

xnode = [
    -8 0;
    -4 2;
    0 4;
    0 0;
];

icone = [1 2;
        1 4;
        2 4;
        2 3;
        3 4;
];

Fixnodes = [
    1 2 0;
    3 1 0;
    4 1 0;
];

Pointload = [
    2 0 -20e3;
    3 0 -10e3;
];

model.nnodes = size(xnode,1);
model.nbarras = size(icone,1);

propiedades(1) = struct('E', 210e9, 'A', 10e-4);
propiedades(2) = struct('E', 210e9, 'A', 10e-4);
propiedades(3) = struct('E', 210e9, 'A', 10e-4);
propiedades(4) = struct('E', 210e9, 'A', 10e-4);
propiedades(5) = struct('E', 210e9, 'A', 5e-4);

[K,F] = barras_initialize(model.nnodes);
[K,F] = barras_gen_system(K,F,xnode,icone,propiedades,model);
[F] = barras_pointload(F,Pointload);
[Kd,Fd] = barras_fixnodes(K,F,Fixnodes);
U = Kd\Fd;
[Def,Ten] = barras_DT(xnode,icone,propiedades, model,U);
reaction = K*U;

function graficarArmadura(xnode, icone, nbarras, U, nodo)
   hold on
   
    for barra = 1 : nbarras
        indiceNodo1 = icone(barra, 1);
        indiceNodo2 = icone(barra, 2);
        nodo1 = xnode(indiceNodo1, :);
        nodo2 = xnode(indiceNodo2, :);
        x = [nodo1(1) nodo2(1)]';
        y = [nodo1(2) nodo2(2)]';
        plot(x, y, "bo-")

        indiceGlobalU1 = indiceNodo1*2 - 1;
        indiceGlobalV1 = indiceNodo1*2;
        indiceGlobalU2 = indiceNodo2*2 - 1;
        indiceGlobalV2 = indiceNodo2*2;
        xNuevo = x + U([indiceGlobalU1 indiceGlobalU2]);
        yNuevo = y + U([indiceGlobalV1 indiceGlobalV2]);
        plot(xNuevo, yNuevo, "ro-")
    end
    
    hold off
    
    title(["Nodo " num2str(nodo)])
    xlabel("x [metros]")
    ylabel("y [metros]")
    grid on
    
    xNodo = xnode(nodo, 1);
    yNodo = xnode(nodo, 2);
    deltaLim = 1.5e-2;
    xlim([xNodo-deltaLim xNodo+deltaLim])
    ylim([yNodo-deltaLim yNodo+deltaLim])
end

puntosDerecha = xnode([2 1], :);
puntosDerecha(:, 1) = -puntosDerecha(:, 1);
uPuntosDerecha = U([3 4 1 2]);
uPuntosDerecha([1 3], 1) = -uPuntosDerecha([1 3], 1);
xnodeAmpliado = [xnode; puntosDerecha];
uAmpliado = [U; uPuntosDerecha];
iconeAmpliado = [
    icone;
    3 5;
    4 5;
    4 6;
    5 6;
];
modelAmpliado.nnodes = size(xnodeAmpliado, 1);
modelAmpliado.nbarras = size(iconeAmpliado, 1);

figure
for nodo = 1 : modelAmpliado.nnodes
    subplot(2, 3, nodo)
    graficarArmadura(xnodeAmpliado, iconeAmpliado, modelAmpliado.nbarras, uAmpliado, nodo)
end

hleg = legend("Configuración inicial", "Configuración final");
newPosition = [0.875 0.8 0.085 0.1];
newUnits = 'normalized';
set(hleg, 'Position', newPosition, 'Units', newUnits);
