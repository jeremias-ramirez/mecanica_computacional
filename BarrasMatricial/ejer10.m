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
    subplot(3, 2, nodo)
    hold on
   
    for barra = 1 : modelAmpliado.nbarras
        indiceNodo1 = iconeAmpliado(barra, 1);
        indiceNodo2 = iconeAmpliado(barra, 2);
        nodo1 = xnodeAmpliado(indiceNodo1, :);
        nodo2 = xnodeAmpliado(indiceNodo2, :);
        x = [nodo1(1) nodo2(1)]';
        y = [nodo1(2) nodo2(2)]';
        plot(x, y, "bo-");

        indiceGlobalU1 = indiceNodo1*2 - 1;
        indiceGlobalV1 = indiceNodo1*2;
        indiceGlobalU2 = indiceNodo2*2 - 1;
        indiceGlobalV2 = indiceNodo2*2;
        xNuevo = x + uAmpliado([indiceGlobalU1 indiceGlobalU2]);
        yNuevo = y + uAmpliado([indiceGlobalV1 indiceGlobalV2]);
        plot(xNuevo, yNuevo, "ro-");
    end
    
    xNodo = xnodeAmpliado(nodo, 1);
    yNodo = xnodeAmpliado(nodo, 2);
    quiver(xNodo, yNodo, uAmpliado(nodo*2-1), uAmpliado(nodo*2), "k", "linewidth", 2, "maxheadsize", 0.1);
    
    hold off
    
    title(["Nodo " num2str(nodo)])
    xlabel("x [metros]")
    ylabel("y [metros]")
    deltaLim = 1.5e-2;
    xlim([xNodo-deltaLim xNodo+deltaLim])
    ylim([yNodo-deltaLim yNodo+deltaLim])
    grid on
end

hleg = legend("Configuración inicial", "Configuración final");
newPosition = [0.82 0.85 0.16 0.05];
newUnits = 'normalized';
set(hleg, 'Position', newPosition, 'Units', newUnits);
%print -dsvg "-S700, 1000" ejer10.svg