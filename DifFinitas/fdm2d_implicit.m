function [PHI_vec, Q_vec] = fdm2d_implicit(K,F,xnode,neighb,model,dt)
    a = (model.rho*model.cp) / dt;
    I = eye(model.nnodes,model.nnodes);
    
    PHI = model.PHI_n; %solucion inicial
    PHI_n = model.PHI_n;
    PHI_vec = PHI;
    Q_vec = zeros(model.nnodes,2);

    for n = 1 : model.maxit
     	
	KK = a*I + K;
        FF = F + a*PHI; 
        PHI = KK\FF;
        
   
        % Error relativo entre las últimas dos iteraciones
        err = norm(PHI-PHI_n,2)/norm(PHI,2);
	
	PHI_n = PHI;
        PHI_vec = [PHI_vec PHI];
        [Q] = fdm2d_flux(PHI,neighb,xnode,model.k);
        Q_vec = [Q_vec Q];

        
        if err < model.tol
          disp('Método terminado por tolerancia de error.');
          return;
        end
    end

    disp('Método terminado por límite de iteraciones.');
end
