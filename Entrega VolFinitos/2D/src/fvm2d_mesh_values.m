function [Z] = fvm2d_mesh_values(Z,X,Y,xnode,icone,neighb,cells,DIR)
    NN = size(icone,1);
    
    Z = full(Z);

    % For each cell centroid we draw triangles to connect boundaries
    for i = 1 : NN
        X(i) = cells(i).cx;
        Y(i) = cells(i).cy;

        if (~isnan(cells(i).ts))
            Z(end+1) =  cells(i).ts;
            E = neighb(i,2);

            if E ~= -1 && ~isnan(cells(E).ts)
                Z(end+1) =  cells(E).ts;
            else
                tc = (cells(i).ts + cells(i).te)/2;
                Z(end+1) = tc;
                Z(end+1) = cells(i).te;
            end
        end

        if (~isnan(cells(i).te))
            Z(end+1) =  cells(i).te;
            N = neighb(i,3);

            if N ~= -1 && ~isnan(cells(N).te)
                Z(end+1) =  cells(N).te;
            else
                tc = (cells(i).te + cells(i).tn)/2;
                Z(end+1) = tc;
                Z(end+1) = cells(i).tn;
            end
        end

        if (~isnan(cells(i).tn))
            Z(end+1) =  cells(i).tn;
            W = neighb(i,4);

            if W ~= -1 && ~isnan(cells(W).tn)
                Z(end+1) =  cells(W).tn;
            else
                tc = (cells(i).tn + cells(i).tw)/2;
                Z(end+1) = tc;
                Z(end+1) = cells(i).tw;
            end
        end

        if (~isnan(cells(i).tw))
            Z(end+1) =  cells(i).tw;
            S = neighb(i,1);

            if S ~= -1 && ~isnan(cells(S).tw)
                Z(end+1) =  cells(S).tw;
            else
                tc = (cells(i).tw + cells(i).ts)/2;
                Z(end+1) = tc;
                Z(end+1) = cells(i).ts;
            end
        end

        % Final step: patches on inner corners (specially for holes)
        E = neighb(i,2);
        N = neighb(i,3);

        if E ~= -1 && N ~= -1 && ~isnan(cells(E).tn) && ~isnan(cells(N).te)
           Z(end+1) = cells(E).tn;
           Z(end+1) = cells(N).te;
           Z(end+1) = (cells(N).te + cells(E).tn)/2;
        end

        W = neighb(i,4);
        N = neighb(i,3);

        if W ~= -1 && N ~= -1 && ~isnan(cells(W).tn) && ~isnan(cells(N).tw)
           Z(end+1) = cells(N).tw;
           Z(end+1) = cells(W).tn;
           Z(end+1) = (cells(N).tw + cells(W).tn)/2;
        end

        W = neighb(i,4);
        S = neighb(i,1);

        if W ~= -1 && S ~= -1 && ~isnan(cells(W).ts) && ~isnan(cells(S).tw)
           Z(end+1) = cells(W).ts;
           Z(end+1) = cells(S).tw;
           Z(end+1) = (cells(S).tw + cells(W).ts)/2;
        end

        E = neighb(i,2);
        S = neighb(i,1);

        if E ~= -1 && S ~= -1 && ~isnan(cells(E).ts) && ~isnan(cells(S).te)
           Z(end+1) = cells(E).ts;
           Z(end+1) = cells(S).te;
           Z(end+1) = (cells(S).te + cells(E).ts)/2;
        end
    end
    
    % Final patch for every Dirichlet corner...gives a better visual look
    for i = 1 : size(DIR,1)
        P = DIR(i,1);
        val = DIR(i,2);
        face = DIR(i,3);
        
        P1 = icone(P,1);
        P2 = icone(P,2);
        P3 = icone(P,3);
        P4 = icone(P,4);

        if (face == 1)
            idx = X == xnode(P1,1) & Y == xnode(P1,2);
            Z(idx) = val;
            
            idx = X == xnode(P2,1) & Y == xnode(P2,2);
            Z(idx) = val;
        end
        
        if (face == 2)
            idx = X == xnode(P2,1) & Y == xnode(P2,2);
            Z(idx) = val;
            
            idx = X == xnode(P3,1) & Y == xnode(P3,2);
            Z(idx) = val;
        end
        
        if (face == 3)
            idx = X == xnode(P3,1) & Y == xnode(P3,2);
            Z(idx) = val;
            
            idx = X == xnode(P4,1) & Y == xnode(P4,2);
            Z(idx) = val;
        end
        
        if (face == 4)
            idx = X == xnode(P1,1) & Y == xnode(P1,2);
            Z(idx) = val;
            
            idx = X == xnode(P4,1) & Y == xnode(P4,2);
            Z(idx) = val;
        end
    end
end
