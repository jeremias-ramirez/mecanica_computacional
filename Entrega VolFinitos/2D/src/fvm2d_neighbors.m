function [neighb] = fvm2d_neighbors(icone)
    M = size(icone,1);
    neighb = -1*ones(M,4);

    for i = 1 : M
        P1 = icone(i,1);
        P2 = icone(i,2);
        P3 = icone(i,3);
        P4 = icone(i,4);

        S = find(icone(:,4) == P1 & icone(:,3) == P2);
        E = find(icone(:,1) == P2 & icone(:,4) == P3);
        N = find(icone(:,1) == P4 & icone(:,2) == P3);
        W = find(icone(:,2) == P1 & icone(:,3) == P4);
        
        if ~isempty(S)
            neighb(i,1) = S;
        end
        
        if ~isempty(E)
            neighb(i,2) = E;
        end
        
        if ~isempty(N)
            neighb(i,3) = N;
        end
        
        if ~isempty(W)
            neighb(i,4) = W;
        end
    end
end