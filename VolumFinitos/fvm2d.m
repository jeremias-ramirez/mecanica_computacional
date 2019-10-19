function [PHI,Q] = fvm2d(xnode,icone,DIR,NEU,ROB,model)
    % - xnode : matriz de pares (x,y). Cada fila contiene las coordenadas
    %   x e y de un nodo de la malla. El número de fila (índice) corresponde
    %   al número de nodo.
    % - icone : matriz de conectividad. Cada fila contiene los índices de
    %   los nodos que integran un elemento, típicamente los 4 nodos que forman
    %   un cuadrángulo.
    
    % Armado de la matriz de vecindad
    [neighb] = fvm2d_neighbors(icone);
    
    % Inicialización de variables principales del sistema
    [K,F,cells] = fvm2d_initialize(xnode,icone,neighb,model.th,model.k);
    
    % Ensamble de coeficientes del sistema
    [K,F] = fvm2d_gen_system(K,F,neighb,cells,model.c,model.G);

    % Ensamble de celdas con lados en frontera Neumann
    [F] = fvm2d_neumann(F,cells,NEU);
    
    % Ensamble de celdas con lados en  frontera Robin
    [K,F] = fvm2d_robin(K,F,cells,ROB);    
    
    % Ensamble de celdas con lados en  frontera Dirichlet
    [K,F] = fvm2d_dirichlet(K,F,cells,DIR);
    
    % Resolución del sistema lineal de ecuaciones
    [PHI,Q] = fvm2d_solve(K,F,neighb,cells,model);
end

