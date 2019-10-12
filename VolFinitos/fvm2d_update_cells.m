function [cells] = fvm2d_update_cells(cells,DIR,NEU,ROB,PHI)
    %% Definition of function to calculate Robin coefficients
    function [a,b,c,d] = fvm2d_robin_coeff(xP,xB,k,h,phi_inf)
        delta = xP - xB;
        a = (h*phi_inf)/(k - h*delta);
        b = (-h)/(k - h*delta);

        c = -delta*a;
        d = 1 - delta*b;
    end

    %% Update Dirichlet borders
    M = size(DIR,1);
    for i = 1 : M
        P = DIR(i,1);
        val = DIR(i,2);
        
        if DIR(i,3) == 1
            cells(P).ts = val;
        end
        
        if DIR(i,3) == 2
            cells(P).te = val;
        end
        
        if DIR(i,3) == 3
            cells(P).tn = val;
        end
        
        if DIR(i,3) == 4
            cells(P).tw = val;
        end
    end
    
    %% Update Neumann borders
    M = size(NEU,1);
    for i = 1 : M
        P = NEU(i,1);
        q = NEU(i,2);
        
        if NEU(i,3) == 1
            cells(P).ts = PHI(P) - q*cells(P).ds/cells(P).ks;
        end
        
        if NEU(i,3) == 2
            cells(P).te = PHI(P) - q*cells(P).de/cells(P).ke;
        end
        
        if NEU(i,3) == 3
            cells(P).tn = PHI(P) - q*cells(P).dn/cells(P).kn;
        end
        
        if NEU(i,3) == 4
            cells(P).tw = PHI(P) - q*cells(P).dw/cells(P).kw;
        end
    end
    
    %% Update Robin borders
    M = size(ROB,1);
    for i = 1 : M
        P = ROB(i,1);
        h = ROB(i,2);
        phi_inf = ROB(i,3);
        
        if ROB(i,4) == 1
            xP = cells(P).cy;
            xB = cells(P).cy + cells(P).dy/2;  %% REVISAR SIGNOS Y NORMALES
            [~,~,a,b] = fvm2d_robin_coeff(xP,xB,cells(P).ks,h,phi_inf);
            cells(P).ts = a + b*PHI(P);
        end
        
        if ROB(i,4) == 2
            xP = cells(P).cx;
            xB = cells(P).cx + cells(P).dx/2;
            [~,~,a,b] = fvm2d_robin_coeff(xP,xB,cells(P).ke,h,phi_inf);
            cells(P).te = a + b*PHI(P);
        end
        
        if ROB(i,4) == 3
            xP = cells(P).cy;
            xB = cells(P).cy + cells(P).dy/2;
            [~,~,a,b] = fvm2d_robin_coeff(xP,xB,cells(P).kn,h,phi_inf);
            cells(P).tn = a + b*PHI(P);
        end
        
        if ROB(i,4) == 4
            xP = cells(P).cx;
            xB = cells(P).cx + cells(P).dx/2; %% REVISAR SIGNOS Y NORMALES
            [~,~,a,b] = fvm2d_robin_coeff(xP,xB,cells(P).kw,h,phi_inf);
            cells(P).tw = a + b*PHI(P);
        end
    end
    
    %% update cells' centroids
    NN = size(PHI,1);   % number of centroids
    for i = 1 : NN
        cells(i).tc = PHI(i);
    end
end

