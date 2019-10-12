function [X,Y,tri,rec] = fvm2d_mesh_coord(xnode,icone,neighb,cells)
    NN = size(icone,1);

    X = zeros(NN,1);
    Y = zeros(NN,1);
    
    tri = [];
    rec = [];

    %% First we find triangles between cell centroids
    for i = 1 : length(cells)
        P = i;
        E = neighb(P,2);
        N = neighb(P,3);

        if E ~= -1
            NE = neighb(E,3);
            if NE ~= -1 && N ~= -1
                tri = [tri; P E NE];%#ok<*AGROW>
                tri = [tri; P NE N]; 
                rec = [rec; P E NE N];
            end
        end
    end

    %% Then, for each cell centroid we draw triangles to connect boundaries
    for i = 1 : NN
        X(i) = cells(i).cx;
        Y(i) = cells(i).cy;

        if (~isnan(cells(i).ts))
            P1 = icone(i,1);

            coord_x = xnode(P1,1) + cells(i).dx/2;
            coord_y = xnode(P1,2);

            X(end+1) = coord_x;
            Y(end+1) = coord_y;

            E = neighb(i,2);

            if E ~= -1 && ~isnan(cells(E).ts)
                N1 = length(X);
                N2 = E;
                N3 = i;

                P1 = icone(E,1);

                coord_x = xnode(P1,1) + cells(E).dx/2;
                coord_y = xnode(P1,2);

                X(end+1) = coord_x;
                Y(end+1) = coord_y;

                N4 = N1+1;

                tri = [tri; N1 N2 N3;];
                tri = [tri; N1 N4 N2;];
                rec = [rec; N1 N4 N2 N3];
            else
                P2 = icone(i,2);

                M = length(X);

                X(end+1) = xnode(P2,1);
                Y(end+1) = xnode(P2,2);

                X(end+1) = xnode(P2,1);
                Y(end+1) = xnode(P2,2) + cells(i).dy/2;

                N1 = M;
                N2 = M + 1;
                N3 = M + 2;
                N4 = i;

                tri = [tri; N2 N3 N4;];
                tri = [tri; N1 N2 N4;];
                rec = [rec; N2 N3 N4 N1];
            end
        end

        if (~isnan(cells(i).te))
            P2 = icone(i,2);

            X(end+1) = xnode(P2,1);
            Y(end+1) = xnode(P2,2) + cells(i).dy/2;

            N = neighb(i,3);

            if N ~= -1 && ~isnan(cells(N).te)
                N1 = i;
                N2 = N;
                N3 = size(X,1);

                P2 = icone(N,2);

                X(end+1) = xnode(P2,1);
                Y(end+1) = xnode(P2,2) + cells(N).dy/2;

                N4 = N3 + 1;

                tri = [tri; N1 N4 N2;];
                tri = [tri; N3 N4 N1;];
                rec = [rec; N1 N2 N4 N3];
            else
                P3 = icone(i,3);

                M = length(X);

                X(end+1) = xnode(P3,1);
                Y(end+1) = xnode(P3,2);

                X(end+1) = xnode(P3,1) - cells(i).dx/2;
                Y(end+1) = xnode(P3,2);

                N1 = M;
                N2 = M + 1;
                N3 = M + 2;
                N4 = i;

                tri = [tri; N2 N4 N3;];
                tri = [tri; N1 N2 N4;];
                rec = [rec; N2 N3 N4 N1];
            end
        end

        if (~isnan(cells(i).tn))
            P4 = icone(i,4);

            coord_x = xnode(P4,1) + cells(i).dx/2;
            coord_y = xnode(P4,2);

            X(end+1) = coord_x;
            Y(end+1) = coord_y;

            W = neighb(i,4);

            if W ~= -1 && ~isnan(cells(W).tn)
                N1 = length(X);
                N2 = W;
                N3 = i;

                P4 = icone(W,4);

                coord_x = xnode(P4,1) + cells(W).dx/2;
                coord_y = xnode(P4,2);

                X(end+1) = coord_x;
                Y(end+1) = coord_y;

                N4 = N1+1;

                tri = [tri; N1 N2 N3;];
                tri = [tri; N1 N4 N2;];
                rec = [rec; N1 N4 N2 N3];
            else
                P4 = icone(i,4);

                M = length(X);

                X(end+1) = xnode(P4,1);
                Y(end+1) = xnode(P4,2);

                X(end+1) = xnode(P4,1);
                Y(end+1) = xnode(P4,2) - cells(i).dy/2;

                N1 = M;
                N2 = M + 1;
                N3 = M + 2;
                N4 = i;

                tri = [tri; N2 N3 N4;];
                tri = [tri; N1 N2 N4;];
                rec = [rec; N1 N2 N3 N4];
            end
        end

        if (~isnan(cells(i).tw))
            P1 = icone(i,1);

            coord_x = xnode(P1,1);
            coord_y = xnode(P1,2) + cells(i).dy/2;

            X(end+1) = coord_x;
            Y(end+1) = coord_y;

            S = neighb(i,1);

            if S ~= -1 && ~isnan(cells(S).tw)
                N1 = length(X);
                N2 = S;
                N3 = i;

                P1 = icone(S,1);

                coord_x = xnode(P1,1);
                coord_y = xnode(P1,2) + cells(S).dy/2;

                X(end+1) = coord_x;
                Y(end+1) = coord_y;

                N4 = N1+1;

                tri = [tri; N2 N3 N4;];
                tri = [tri; N1 N4 N3;];
                rec = [rec; N2 N3 N1 N4];
            else
                tc = (cells(i).tw + cells(i).ts)/2;
                P1 = icone(i,1);

                M = length(X);

                X(end+1) = xnode(P1,1);
                Y(end+1) = xnode(P1,2);

                X(end+1) = xnode(P1,1) + cells(i).dx/2;
                Y(end+1) = xnode(P1,2);

                N1 = M;
                N2 = M + 1;
                N3 = M + 2;
                N4 = i;

                tri = [tri; N2 N3 N4;];
                tri = [tri; N1 N2 N4;];
                rec = [rec; N2 N3 N4 N1];
            end
        end

        %% Final step: patches on inner corners (specially for holes)
        E = neighb(i,2);
        N = neighb(i,3);

        if E ~= -1 && N ~= -1 && ~isnan(cells(E).tn) && ~isnan(cells(N).te)
           Nc = icone(i,3);

           N1 = i;
           N2 = E;
           N3 = N;

           X(end+1) = xnode(Nc,1) + cells(E).dx/2;
           Y(end+1) = xnode(Nc,2);

           N4 = length(X);

           X(end+1) = xnode(Nc,1);
           Y(end+1) = xnode(Nc,2) + cells(N).dy/2;

           N5 = N4 + 1;

           X(end+1) = xnode(Nc,1);
           Y(end+1) = xnode(Nc,2);

           N6 = N5 + 1;

           tri = [tri; N1 N2 N6; N2 N4 N6];
           tri = [tri; N1 N6 N3; N6 N5 N3];
           rec = [rec; N1 N2 N4 N6];
           rec = [rec; N3 N1 N6 N5];
        end

        W = neighb(i,4);
        N = neighb(i,3);

        if W ~= -1 && N ~= -1 && ~isnan(cells(W).tn) && ~isnan(cells(N).tw)
           Nc = icone(i,4);

           N1 = i;
           N2 = N;
           N3 = W;

           X(end+1) = xnode(Nc,1);
           Y(end+1) = xnode(Nc,2) + cells(N).dy/2;

           N4 = length(X);

           X(end+1) = xnode(Nc,1) - cells(W).dx/2;
           Y(end+1) = xnode(Nc,2);

           N5 = N4 + 1;

           X(end+1) = xnode(Nc,1);
           Y(end+1) = xnode(Nc,2);

           N6 = N5 + 1;

           tri = [tri; N1 N2 N6; N2 N4 N6];
           tri = [tri; N1 N6 N3; N6 N5 N3];
           rec = [rec; N1 N2 N4 N6];
           rec = [rec; N3 N1 N6 N5];
        end

        W = neighb(i,4);
        S = neighb(i,1);

        if W ~= -1 && S ~= -1 && ~isnan(cells(W).ts) && ~isnan(cells(S).tw)
           Nc = icone(i,1);

           N1 = i;
           N2 = W;
           N3 = S;

           X(end+1) = xnode(Nc,1) - cells(W).dx/2;
           Y(end+1) = xnode(Nc,2);

           N4 = length(X);

           X(end+1) = xnode(Nc,1);
           Y(end+1) = xnode(Nc,2) - cells(S).dy/2;

           N5 = N4 + 1;

           X(end+1) = xnode(Nc,1);
           Y(end+1) = xnode(Nc,2);

           N6 = N5 + 1;

           tri = [tri; N1 N2 N6; N2 N4 N6];
           tri = [tri; N1 N6 N3; N6 N5 N3];
           rec = [rec; N1 N2 N4 N6];
           rec = [rec; N3 N1 N6 N5];
        end

        E = neighb(i,2);
        S = neighb(i,1);

        if E ~= -1 && S ~= -1 && ~isnan(cells(E).ts) && ~isnan(cells(S).te)
           Nc = icone(i,2);

           N1 = i;
           N2 = E;
           N3 = S;

           X(end+1) = xnode(Nc,1) + cells(E).dx/2;
           Y(end+1) = xnode(Nc,2);

           N4 = length(X);

           X(end+1) = xnode(Nc,1);
           Y(end+1) = xnode(Nc,2) - cells(S).dy/2;

           N5 = N4 + 1;

           X(end+1) = xnode(Nc,1);
           Y(end+1) = xnode(Nc,2);

           N6 = N5 + 1;

           tri = [tri; N1 N2 N6; N2 N4 N6];
           tri = [tri; N1 N6 N3; N6 N5 N3];
           rec = [rec; N1 N2 N4 N6];
           rec = [rec; N3 N1 N6 N5];
        end
    end
end

