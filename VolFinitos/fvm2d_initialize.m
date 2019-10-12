function [K,F,cells] = fvm2d_initialize(xnode,icone,neighb,th,k)
    %% Cantidad de celdas
    NN = size(icone,1);

    %% Inicialización de la matriz K y el vector F (RHS)
    K = sparse(NN);
    F = sparse(NN,1);


    %% Definición de celdas. Arreglo de NN structs para almacenar los parámetros 
    %  de cada celda
    cells = struct('ds', {}, 'de', {}, 'dn', {}, 'dw', {},...
                   'as', {}, 'ae', {}, 'an', {}, 'aw', {},...
                   'ks', {}, 'ke', {}, 'kn', {}, 'kw', {},...
                   'ts', {}, 'te', {}, 'tn', {}, 'tw', {}, 'tc', {},...
                   'dx', {}, 'dy', {}, 'v',{}, 'cx', {}, 'cy', {});
    
    %% Bucle principal para especificar todos los valores de las celdas
    for i = 1 : NN
        % Vecinos de la celda (i)
        S = neighb(i,1);
        E = neighb(i,2);
        N = neighb(i,3);
        W = neighb(i,4);
        
        % Nodos relevantes para obtener la geometría de la celda (i)
        P1 = icone(i,1);
        P2 = icone(i,2);
        P3 = icone(i,3);
        
        %% Inicializar dx (ancho), dy (alto) y v (volumen) de la celda
        dx = abs(xnode(P2,1) - xnode(P1,1));
        dy = abs(xnode(P3,2) - xnode(P2,2));
        v = dx*dy*th;
        
        cells(i).dx = dx; cells(i).dy = dy; cells(i).v = v;
        
        % Inicializar el centroide de la celda
        cells(i).cx = xnode(P1,1) + dx/2;
        cells(i).cy = xnode(P1,2) + dy/2;
        
        %% Inicializar las distancias entre centroides de celdas
        if S ~= -1
            S2 = icone(S, 2);
            S3 = icone(S, 3);
            dys = abs(xnode(S3,2) - xnode(S2,2));

            ds = dy/2 + dys/2;
        else
            ds = dy/2;
        end

        if E ~= -1
            E1 = icone(E, 1);
            E2 = icone(E, 2);
            dxe = abs(xnode(E2,1) - xnode(E1,1));

            de = dx/2 + dxe/2;
        else
            de = dx/2;
        end

        if N ~= -1
            N2 = icone(N, 2);
            N3 = icone(N, 3);
            dyn = abs(xnode(N3,2) - xnode(N2,2));

            dn = dy/2 + dyn/2;
        else
            dn = dy/2;
        end

        if W ~= -1
            W1 = icone(W, 1);
            W2 = icone(W, 2);
            dxw = abs(xnode(W2,1) - xnode(W1,1));

            dw = dx/2 + dxw/2;
        else
            dw = dx/2;
        end
        
        cells(i).ds = ds; cells(i).de = de; cells(i).dn = dn; cells(i).dw = dw;
        
        %% Inicializar las áreas de las caras de la celda (i)
        as = dx*th;
        ae = dy*th;
        an = dx*th;
        aw = dy*th;
        
        cells(i).as = as; cells(i).ae = ae; cells(i).an = an; cells(i).aw = aw;
        
        %% Inicialización de las conductividades en las caras de la celda (i)
        ks = k(i); ke = k(i); kn = k(i); kw = k(i);

        if(S ~= -1)
            if (k(i) ~= k(S))
                S2 = icone(S, 2);
                S3 = icone(S, 3);
                dS = abs(xnode(S3,2) - xnode(S2,2))/2;
                dP = dy/2;
                f = dS/(dP + dS);

                ks = k(P)*k(S)/(f*k(P) + (1-f)*k(S));
            end
        end

        if(E ~= -1)
            if (k(i) ~= k(E))
                E1 = icone(E, 1);
                E2 = icone(E, 2);
                dE = abs(xnode(E2,1) - xnode(E1,1))/2;
                dP = dx/2;
                f = dE/(dP + dE);

                ke = k(P)*k(E)/(f*k(P) + (1-f)*k(E));
            end
        end

        if(N ~= -1)
            if (k(i) ~= k(N))
                N2 = icone(N, 2);
                N3 = icone(N, 3);
                dN = abs(xnode(N3,2) - xnode(N2,2))/2;
                dP = dy/2;
                f = dN/(dP + dN);

                kn = k(P)*k(N)/(f*k(P) + (1-f)*k(N));
            end
        end

        if(W ~= -1)
            if (k(i) ~= k(W))
                W1 = icone(W, 1);
                W2 = icone(W, 2);
                dW = abs(xnode(W2,1) - xnode(W1,1))/2;
                dP = dx/2;
                f = dW/(dP + dW);

                kw = k(P)*k(W)/(f*k(P) + (1-f)*k(W));
            end
        end
        
        cells(i).ks = ks; cells(i).ke = ke; cells(i).kn = kn; cells(i).kw = kw;
        
        %% Inicialización de las temperaturas de la celda como NaN (útil para graficar)
        cells(i).ts = nan; cells(i).te = nan; % valores al sur y al este, respectivamente
        cells(i).tn = nan; cells(i).tw = nan; % valores al norte y al oeste, respectivamente
        cells(i).tc = nan; % valor en el centroide
    end
end
