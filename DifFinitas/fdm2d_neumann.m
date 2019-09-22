function [F] = fdm2d_neumann(F,xnode,neighb,NEU)
    M = size(NEU, 1);
    for n = 1 : M
        P = NEU(n, 1);           
        S = neighb(P, 1);        
        E = neighb(P, 2);       
        N = neighb(P, 3);        
        W = neighb(P, 4);        

        q = NEU(n,2);
                
        if (NEU(n,3) == 1)
            dy = abs(xnode(N,2) - xnode(P,2));
            F(P) = F(P) - 2*q/dy;
        end
        
        if (NEU(n,3) == 2) 
            dx = abs(xnode(W,1) - xnode(P,1));
            F(P) = F(P) - 2*q/dx;
        end
        
        if (NEU(n,3) == 3) 
            dy = abs(xnode(S,2) - xnode(P,2));
            F(P) = F(P) - 2*q/dy;
        end
        
        if (NEU(n,3) == 4)
            dx = abs(xnode(E,1) - xnode(P,1));
            F(P) = F(P) - 2*q/dx;
        end
    end
end