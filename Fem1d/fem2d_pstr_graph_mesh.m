function [] = fem2d_pstr_graph_mesh(U,reaction,Ten_VM,xnode,icone,mode,graph)
    figure('Name', 'Resultados');
    
    if (graph < 0 || graph > 5)
        graph = 0;
    end
    
    triangles = [];
    rectangles = [];
    for i = 1 : size(icone,1)
        if icone(i,4) == -1
            triangles = [triangles; icone(i,1:3)]; %#ok<*AGROW>
        else
            rectangles = [rectangles; icone(i,:)];
        end
    end

    trimesh = triangles;
    for i = 1 : size(icone,1)
        if icone(i,4) ~= -1
            i1 = icone(i,1);
            i3 = icone(i,3);
            i4 = icone(i,4);
            trimesh = [trimesh; icone(i,1:3); i1,i3,i4];
            rectangles = [rectangles; icone(i,:)];
        end
    end

    %% Displacements
    if graph == 0
        clf;
        hold on;
        axis equal;

        Z = zeros(size(xnode,1),1);

        patch('Vertices',xnode,'Faces',triangles,'FaceVertexCData',Z,'FaceColor','none','EdgeColor','red');
        patch('Vertices',xnode,'Faces',rectangles,'FaceVertexCData',Z,'FaceColor','none','EdgeColor','red');

        xnode_d = xnode;
        xnode_d(:,1) = xnode_d(:,1) + U(1:2:length(U));
        xnode_d(:,2) = xnode_d(:,2) + U(2:2:length(U));

        quiver(xnode(:,1),xnode(:,2),U(1:2:length(U)),U(2:2:length(U)),'color', 'black');

        patch('Vertices',xnode_d,'Faces',triangles,'FaceVertexCData',Z,'FaceColor','none','EdgeColor','black');
        patch('Vertices',xnode_d,'Faces',rectangles,'FaceVertexCData',Z,'FaceColor','none','EdgeColor','black');
        hold off;
    end

    %% Von Misses (scalar)
    if graph == 1
        X = xnode(:,1);
        Y = xnode(:,2);
        Z = Ten_VM;

        xmax = max(X);
        xmin = min(X);
        ymax = max(Y);
        ymin = min(Y);
        zmax = max(Z);
        zmin = min(Z);

        ratio = (xmax-xmin)/(ymax-ymin);

        clf;
        	
        hold on;
        patch('Vertices',[X Y Z],'Faces',trimesh,'FaceVertexCData',Z);
        hold off;
        shading interp;

        if (mode == 0)
            hold on;
            patch('Vertices',[X Y Z],'Faces',triangles,'FaceVertexCData',Z,...
            'FaceColor','none');

            patch('Vertices',[X Y Z],'Faces',rectangles,'FaceVertexCData',Z,...
            'FaceColor','none');
            hold off;
        end

        view(2);
        grid on;
        pbaspect([ratio 1 1]);
        zlim([zmin-0.001, zmax+0.001]);
        drawnow;            
        colorbar;
    end

    %% Reaction (vectorial)
    if graph == 2
        clf;
        hold on;
        axis equal;

        Z = zeros(size(xnode,1),1);

        patch('Vertices',xnode,'Faces',triangles,'FaceVertexCData',Z,'FaceColor','none','EdgeColor','black');
        patch('Vertices',xnode,'Faces',rectangles,'FaceVertexCData',Z,'FaceColor','none','EdgeColor','black');

        quiver(xnode(:,1),xnode(:,2),reaction(1:2:length(reaction)),reaction(2:2:length(reaction)),'color', 'black');
        hold off;
    end
end

