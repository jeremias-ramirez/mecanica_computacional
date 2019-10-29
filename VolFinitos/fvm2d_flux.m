function [Q] = fvm2d_flux(PHI,cells,neighb)
    Qx = zeros(size(neighb,1),1);
    Qy = zeros(size(neighb,1),1);
    
    PHI = full(PHI);
    
    cx = zeros(size(neighb,1),1);
    cy = zeros(size(neighb,1),1);

    for P = 1 : size(neighb,1)
        S = neighb(P,1);
        E = neighb(P,2);
        N = neighb(P,3);
        W = neighb(P,4);

        if (E ~= -1 && W ~= -1)
            dE = cells(P).de;
            dW = cells(P).dw;
            fx = dE / (dE + dW);
            k = fx*cells(P).ke + (1-fx)*cells(P).kw;
            Qx(P) = -k * (PHI(E) - PHI(W)) / (dE + dW);
        elseif (E ~= -1)
            Qx(P) = -cells(P).ke * (PHI(E) - PHI(P)) / cells(P).de;
        else
            Qx(P) = -cells(P).kw * (PHI(P) - PHI(W)) / cells(P).dw;
        end
        
        if (N ~= -1 && S ~= -1)
            dN = cells(P).dn;
            dS = cells(P).ds;
            fx = dN / (dN + dS);
            k = fx*cells(P).ks + (1-fx)*cells(P).kn;
            Qy(P) = -k * (PHI(N) - PHI(S)) / (dN + dS);
        elseif (N ~= -1)
            Qy(P) = -cells(P).kn * (PHI(N) - PHI(P)) / cells(P).dn;
        else
            Qy(P) = -cells(P).ks * (PHI(P) - PHI(S)) / cells(P).ds;
        end

        cx(P) = cells(P).cx;
        cy(P) = cells(P).cy;
    end
    
    Q = [Qx Qy];
end

