function [PHI,Q] = fvm2d_solve(K,F,neighb,cells,model)
    % Time Scheme: [0] Explicit, [1] Implicit, [X] Stationary
    if model.ts == 0 % Explicit
        disp('Iniciando esquemas temporales...');
        % Explicit (Forward Euler) method's time step
        [dt] = fvm2d_explicit_delta_t(cells,model.k,model.rho,model.cp);
        [PHI,Q] = fvm2d_explicit(K,F,cells,neighb,model,dt);
    elseif model.ts == 1 % Implicit
        disp('Iniciando esquemas temporales...');
        % Arbitrary time step
        dt = model.dt;
        [PHI,Q] = fvm2d_implicit(K,F,cells,neighb,model,dt);
    else % Stationary
        disp('Resolución del sistema de ecuaciones...');
        % Linear equations system solve
        PHI = K\F;
        
        % Calculate system heat flux
        [Q] = fvm2d_flux(PHI,cells,neighb);
    end
end