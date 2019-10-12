function [X,Y,Z] = fvm2d_mesh_interp(xnode,icone,neighb,cells,DIR,NEU,ROB)
    %% Triangular shape functions
    function [Ni, Nj, Nk] = aux_ntriang(i,j,k,p)
        xi = i(1); yi = i(2);   % extracts coordinates [xi,yi] of point (i)
        xj = j(1); yj = j(2);   % extracts coordinates [xj,yj] of point (j)
        xk = k(1); yk = k(2);   % extracts coordinates [xk,yk] of point (k)
        xp = p(1); yp = p(2);   % extracts coordinates [xp,yp] of point (p)

        A = det([1 xi yi; 1 xj yj; 1 xk yk])/2; % element's area
        % shape function's coefficients
        ai = (xj*yk - xk*yj); bi = (yj - yk); ci = (xk - xj);
        aj = (xk*yi - xi*yk); bj = (yk - yi); cj = (xi - xk);
        ak = (xi*yj - xj*yi); bk = (yi - yj); ck = (xj - xi);


        Ni = (ai + bi*xp + ci*yp)/(2*A);
        Nj = (aj + bj*xp + cj*yp)/(2*A);
        Nk = (ak + bk*xp + ck*yp)/(2*A);
    end
    
    %% Quadrangular shape functions 
    function [Ni,Nj,Nk,Nl] = aux_nquad(i,j,k,l,p)
        xi = i(1); yi = i(2);   % extracts coordinates [xi,yi] of point (i)
        xj = j(1); yj = j(2);   % extracts coordinates [xj,yj] of point (j)
        xk = k(1); yk = k(2);   % extracts coordinates [xk,yk] of point (k)
        xl = l(1); yl = l(2);   % extracts coordinates [xl,yl] of point (l)
        xp = p(1); yp = p(2);   % extracts coordinates [xp,yp] of point (p)

        s = (xp - xi) / (xj - xi);
        t = (yp - yj) / (yk - yj);

        Ni = (1 - s) * (1 - t);
        Nj = s*(1-t);
        Nk = s*t;
        Nl = (1-s)*t;
    end

    %% Interpolation
	X = xnode(:,1);
	Y = xnode(:,2);
	Z = NaN*zeros(size(xnode,1),1);
            
    for i = 1 : size(icone,1)
		S = neighb(i,1);
		E = neighb(i,2);
		N = neighb(i,3);
		W = neighb(i,4);

        SE = -1; SW = -1;
        if (S ~= -1)
            SW = neighb(S,4);
            SE = neighb(S,2);
        end
        
        NE = -1; NW = -1;
        if (N ~= -1)
            NE = neighb(N,2);
            NW = neighb(N,4);
        end
        
        P1 = icone(i,1);
		P2 = icone(i,2);
		P3 = icone(i,3);
		P4 = icone(i,4);
        
        if (SW ~= -1)
            if (W ~= -1)
                a = [cells(SW).cx, cells(SW).cy];
                b = [cells(S).cx, cells(S).cy];
                c = [cells(i).cx, cells(i).cy];
                d = [cells(W).cx, cells(W).cy];
                
                [Na,Nb,Nc,Nd] = aux_nquad(a,b,c,d,xnode(P1,:));
                Z(P1) = Na*cells(SW).tc + Nb*cells(S).tc + Nc*cells(i).tc + Nd*cells(W).tc;
            else
                a = [cells(i).cx, cells(i).cy];
                b = [cells(SW).cx, cells(SW).cy];
                c = [cells(S).cx, cells(S).cy];
                
                [Na,Nb,Nc] = aux_ntriang(a,b,c,xnode(P1,:));
                
                Z(P1) = Na*cells(i).tc + Nb*cells(SW).tc + Nc*cells(S).tc;
            end
        end
        
        if (SE ~= -1 && isnan(Z(P2)))
            if (E ~= -1)
                a = [cells(S).cx, cells(S).cy];
                b = [cells(SE).cx, cells(SE).cy];
                c = [cells(E).cx, cells(E).cy];
                d = [cells(i).cx, cells(i).cy];
                
                [Na,Nb,Nc,Nd] = aux_nquad(a,b,c,d,xnode(P2,:));
                Z(P2) = Na*cells(S).tc + Nb*cells(SE).tc + Nc*cells(E).tc + Nd*cells(i).tc;
            else
                a = [cells(i).cx, cells(i).cy];
                b = [cells(S).cx, cells(S).cy];
                c = [cells(SE).cx, cells(SE).cy];
                
                [Na,Nb,Nc] = aux_ntriang(a,b,c,xnode(P2,:));
                
                Z(P2) = Na*cells(i).tc + Nb*cells(S).tc + Nc*cells(SE).tc;
            end
        end
        
        if (NE ~= -1 && isnan(Z(P3)))
            if (E ~= -1)
                a = [cells(i).cx, cells(i).cy];
                b = [cells(E).cx, cells(E).cy];
                c = [cells(NE).cx, cells(NE).cy];
                d = [cells(N).cx, cells(N).cy];
                
                [Na,Nb,Nc,Nd] = aux_nquad(a,b,c,d,xnode(P3,:));
                Z(P3) = Na*cells(i).tc + Nb*cells(E).tc + Nc*cells(NE).tc + Nd*cells(N).tc;
            else
                a = [cells(i).cx, cells(i).cy];
                b = [cells(NE).cx, cells(NE).cy];
                c = [cells(N).cx, cells(N).cy];
                
                [Na,Nb,Nc] = aux_ntriang(a,b,c,xnode(P3,:));
                
                Z(P3) = Na*cells(i).tc + Nb*cells(NE).tc + Nc*cells(N).tc;
            end
        end
        
        if (NW ~= -1 && isnan(Z(P4)))
            if (W ~= -1)
                a = [cells(W).cx, cells(W).cy];
                b = [cells(i).cx, cells(i).cy];
                c = [cells(N).cx, cells(N).cy];
                d = [cells(NW).cx, cells(NW).cy];
                
                [Na,Nb,Nc,Nd] = aux_nquad(a,b,c,d,xnode(P4,:));
                Z(P4) = Na*cells(W).tc + Nb*cells(i).tc + Nc*cells(N).tc + Nd*cells(NW).tc;
            else
                a = [cells(i).cx, cells(i).cy];
                b = [cells(N).cx, cells(N).cy];
                c = [cells(NW).cx, cells(NW).cy];
                
                [Na,Nb,Nc] = aux_ntriang(a,b,c,xnode(P4,:));
                
                Z(P4) = Na*cells(i).tc + Nb*cells(N).tc + Nc*cells(NW).tc;
            end
        end
    end
    
    boundary = [];
    
    if (~isempty(NEU))
        boundary = [boundary; NEU(:,1), NEU(:,3)];
    end
    
    if (~isempty(ROB))
        boundary = [boundary; ROB(:,1), ROB(:,4)];
    end
    
    for i = 1 : size(boundary,1)
        P = boundary(i,1);
        face = boundary(i,2);
        
        S = neighb(P,1);
        E = neighb(P,2);
        N = neighb(P,3);
        W = neighb(P,4);
        
        P1 = icone(P,1);
        P2 = icone(P,2);
        P3 = icone(P,3);
        P4 = icone(P,4);
        
        if (face == 1)
            if (W == -1)
                Z(P1) = (cells(P).ts + cells(P).tw)/2;
            elseif (~isnan(cells(W).ts))
                Z(P1) = (cells(P).ts + cells(W).ts)/2;
            end
            
            if (E == -1)
                Z(P2) = (cells(P).ts + cells(P).te)/2;
            elseif (~isnan(cells(E).ts))
                Z(P2) = (cells(P).ts + cells(E).ts)/2;
            end
        end
        
        if (face == 2)
            if (S == -1)
                Z(P2) = (cells(P).ts + cells(P).te)/2;
            elseif (~isnan(cells(S).te))
                Z(P2) = (cells(P).te + cells(S).te)/2;
            end
            
            if (N == -1)
                Z(P3) = (cells(P).tn + cells(P).te)/2;
            elseif (~isnan(cells(N).te))
                Z(P3) = (cells(P).te + cells(N).te)/2;
            end
        end
        
        if (face == 3)
            if (W == -1)
                Z(P4) = (cells(P).tn + cells(P).tw)/2;
            elseif (~isnan(cells(W).tn))
                Z(P4) = (cells(P).tn + cells(W).tn)/2;
            end
            
            if (E == -1)
                Z(P3) = (cells(P).tn + cells(P).te)/2;
            elseif (~isnan(cells(E).tn))
                Z(P3) = (cells(P).tn + cells(E).tn)/2;
            end
        end
        
        if (face == 4)
            if (S == -1)
                Z(P1) = (cells(P).ts + cells(P).tw)/2;
            elseif (~isnan(cells(S).tw))
                Z(P1) = (cells(P).tw + cells(S).tw)/2;
            end
            
            if (N == -1)
                Z(P4) = (cells(P).tn + cells(P).tw)/2;
            elseif (~isnan(cells(N).tw))
                Z(P4) = (cells(P).tw + cells(N).tw)/2;
            end
        end
    end
    
    for i = 1 : size(DIR,1)
        P = DIR(i,1);
        T = DIR(i,2);
        face = DIR(i,3);
        
        P1 = icone(P,1);
        P2 = icone(P,2);
        P3 = icone(P,3);
        P4 = icone(P,4);
        
        if (face == 1)
            Z(P1) = T;
            Z(P2) = T;
        end
        
        if (face == 2)
            Z(P2) = T;
            Z(P3) = T;
        end
        
        if (face == 3)
            Z(P3) = T;
            Z(P4) = T;
        end
        
        if (face == 4)
            Z(P4) = T;
            Z(P1) = T;
        end
    end
end
	