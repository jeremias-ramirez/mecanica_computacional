function [] = fvm2d_graph_mesh(PHI,Q,xnode,icone,neighb,DIR,NEU,ROB,model,mode,graph)
    figure('Name', 'Resultados');
    
    if (mode == 0 || mode == 2 || mode == 4)
        vm = 2;
    elseif (mode == 1 || mode == 3 || mode == 5)
        vm = 3;
    else
        vm = 2;
        mode = 0;
    end
    
    if (graph < 0 || graph > 5)
        graph = 0;
    end
    
    % Initialize cells
    [~,~,cells] = fvm2d_initialize(xnode,icone,neighb,model.th,model.k);
    
    % Update cell's face value (temperature)
    [cells] = fvm2d_update_cells(cells,DIR,NEU,ROB,PHI(:,end));
    
    % Cell's centroids
    cx = zeros(model.ncells,1);
    cy = zeros(model.ncells,1);
    triangles = [];
    rectangles =[];
    for P = 1 : model.ncells
        cx(P) = cells(P).cx;
        cy(P) = cells(P).cy;
        
        E = neighb(P,2);
        N = neighb(P,3);

        if E ~= -1
            NE = neighb(E,3);
            if NE ~= -1 && N ~= -1
                triangles = [triangles; P E NE];%#ok<*AGROW>
                triangles = [triangles; P NE N]; 
                rectangles = [rectangles; P E NE N];
            end
        end
    end
    
    err = 1;

    for i = 1 : size(PHI,2) - 1
        err = [err; norm(PHI(:,i+1)-PHI(:,i),2)/norm(PHI(:,i),2)];
    end

    %% Temperature (scalar)
    if graph == 0
        % Set z-limits for better visualization
        [X,Y,tri,rec] = fvm2d_mesh_coord(xnode,icone,neighb,cells);
        [Z] = fvm2d_mesh_values(PHI(:,end),X,Y,xnode,icone,neighb,cells,DIR);
        zmin = min(Z);
        zmax = max(Z);
        
        for i = 1 : size(PHI,2)
            clf;
            Z = PHI(:,i);
        
            % Update cell's face value (temperature)
            [cells] = fvm2d_update_cells(cells,DIR,NEU,ROB,Z);

            [Z] = fvm2d_mesh_values(Z,X,Y,xnode,icone,neighb,cells,DIR);

            hold on;
            patch('Vertices',[X Y Z],'Faces',tri,'FaceVertexCData',Z);
            hold off;
            shading interp;

            if (mode < 4)
                hold on;
                patch('Vertices',[X Y Z],'Faces',rec,'FaceVertexCData',Z,...
                    'FaceColor','none');
                hold off;
            end

            if (i > 1)
                title(sprintf('nit: %d - error: %e',i-1,err(i)));
            end
            
            if (mode == 2 || mode == 3)
                hold on;
                scatter3(cx,cy,full(PHI(:,i)),200/sqrt(size(PHI(:,i),1)),'k');
                hold off;
            end
            
            view(vm);
            zlim([zmin zmax]);
            grid on;
            drawnow;
            pause(0.000001);
        end
        
        colorbar;
    end
    
    %% Temperature (interpolated)
    if graph == 1
        % Set z-limits for better visualization
        [cells] = fvm2d_update_cells(cells,DIR,NEU,ROB,PHI(:,end));
        [~,~,Z] = fvm2d_mesh_interp(xnode,icone,neighb,cells,DIR,NEU,ROB);
        zmin = min(Z);
        zmax = max(Z);
        
        tri = [];
        for i = 1 : size(icone,1)
            P1 = icone(i,1);
            P2 = icone(i,2);
            P3 = icone(i,3);
            P4 = icone(i,4);
            tri = [tri; P1 P2 P3];
            tri = [tri; P1 P3 P4];
        end
        
        for i = 1 : size(PHI,2)
            clf;
            [cells] = fvm2d_update_cells(cells,DIR,NEU,ROB,PHI(:,i));
            [X,Y,Z] = fvm2d_mesh_interp(xnode,icone,neighb,cells,DIR,NEU,ROB);
            hold on;
            patch('Vertices',[X Y Z],'Faces',tri,'FaceVertexCData',Z);
            hold off;
            shading interp;

            if (mode < 4)
                hold on;
                patch('Vertices',[X Y Z],'Faces',icone,'FaceVertexCData',Z,...
                    'FaceColor','none');
                hold off;
            end

            if (i > 1)
                title(sprintf('nit: %d - error: %e',i-1,err(i)));
            end
            
            if (mode == 2 || mode == 3)
                hold on;
                scatter3(cx,cy,full(PHI(:,i)),200/sqrt(size(PHI(:,i),1)),'k');
                hold off;
            end

            view(vm);
            zlim([zmin*0.999 zmax*1.001]);
            grid on;
            drawnow;
            pause(0.000001);
        end
        
        colorbar;
        
        
    end
    
    %% Heat Flux (vectorial)
    if graph == 2        
        [X,Y,~,rec] = fvm2d_mesh_coord(xnode,icone,neighb,cells);
        [Z] = fvm2d_mesh_values(PHI(:,end),X,Y,xnode,icone,neighb,cells,DIR);
        
        maxq = max(max(abs(Q)));
        
        if (maxq > 0 && maxq < 5)
            scale = 0.5/maxq;
        elseif (maxq >= 5 && maxq < 10)
            scale = 1/maxq;
        elseif (maxq >= 10 && maxq < 100)
            scale = 5/maxq;
        elseif (maxq >= 100 && maxq < 1000)
            scale = 20/maxq;
        elseif (maxq >= 1000)
            scale = 50/maxq;
        else
            scale = 1;
        end
        
        xmax = max(xnode(:,1));
        xmin = min(xnode(:,1));
        ymax = max(xnode(:,2));
        ymin = min(xnode(:,2));
        zmax = max(max(Z));
        zmin = min(min(Z));
        
        for i = 1 : size(PHI,2)
            clf;
            Z = full(PHI(:,i));
        
            % Update cell's face value (temperature)
            [cells] = fvm2d_update_cells(cells,DIR,NEU,ROB,Z);

            [Temp] = fvm2d_mesh_values(Z,X,Y,xnode,icone,neighb,cells,DIR);
            
            if (mode < 4)              
                hold on;
                patch('Vertices',[X Y Temp],'Faces',rec,'FaceVertexCData',Temp,...
                    'FaceColor','none');
                hold off;
                shading interp;
            end 
                
            hold on;
            Qx = Q(:,2*(i-1)+1);
            Qy = Q(:,2*(i-1)+2);
            Qz = zeros(size(Qx,1),1);
            
            quiver3(cx, cy, Z, Qx, Qy, Qz,...
                scale, 'color', 'black');
            hold off;
            if (i > 1)
                title(sprintf('nit: %d - error: %e',i-1,err(i)));
            end

            xlim([xmin-(xmax-xmin)*.1, xmax+(xmax-xmin)*.1]);
            ylim([ymin-(ymax-ymin)*.1, ymax+(ymax-ymin)*.1]);
            zlim([zmin-(zmax-zmin)*.1, zmax+(zmax-zmin)*.1]);
            
            if (mode == 2 || mode == 3)
                hold on;
                scatter3(cx,cy,full(PHI(:,i)),200/sqrt(size(PHI(:,i),1)),'k');
                hold off;
            end
            
            grid on;
            view(vm);
            drawnow;
            pause(0.000001);
        end
        colorbar;
    end
    
    %% Heat Flux x-axis (scalar)
    if graph == 3
        zmax = max(max(Q));
        zmin = min(min(Q));
        
        for i = 1 : size(PHI,2)
            clf;
            Z = Q(:,2*(i-1)+1);
            trisurf(triangles,cx,cy,Z);
            shading interp;

            if (mode < 4)
                hold on;
                patch('Vertices',[cx cy Z],'Faces',rectangles,'FaceVertexCData',Z,...
                'FaceColor','none');
                hold off;
            end
            
            if (i > 1)
                title(sprintf('nit: %d - error: %e',i-1,err(i)));
            end
            
            if (mode == 2 || mode == 3)
                hold on;
                scatter3(cx,cy,Z,200/sqrt(size(Z,1)),'k');
                hold off;
            end

            zlim([zmin zmax]);
            grid on;
            view(vm);
            drawnow;
            pause(0.000001);
        end
        colorbar;
    end
    
    %% Heat Flux y-axis (scalar)
    if graph == 4
        zmax = max(max(Q));
        zmin = min(min(Q));
        
        for i = 1 : size(PHI,2)
            clf;
            Z = Q(:,2*(i-1)+2);
            trisurf(triangles,cx,cy,Z);
            shading interp;

            if (mode < 4)
                hold on;
                patch('Vertices',[cx cy Z],'Faces',rectangles,'FaceVertexCData',Z,...
                'FaceColor','none');
                hold off;
            end
            
            if (i > 1)
                title(sprintf('nit: %d - error: %e',i-1,err(i)));
            end
            
            if (mode == 2 || mode == 3)
                hold on;
                scatter3(cx,cy,Z,200/sqrt(size(Z,1)),'k');
                hold off;
            end

            zlim([zmin zmax]);
            grid on;
            view(vm);
            drawnow;
            pause(0.000001);
        end
        colorbar;
    end
    
    %% Heat Flux Magnitude (scalar)
    if graph == 5
        Z = zeros(size(PHI));
        for i = 1 : size(PHI,2)
            Qx = Q(:,2*(i-1)+1);
            Qy = Q(:,2*(i-1)+2);
            
            for j = 1 : size(Qx,1)
                Z(j,i) = norm([Qx(j), Qy(j)],2);
            end
        end
        
        zmax = max(max(Z));
        zmin = min(min(Z));
        
        for i = 1 : size(PHI,2)
            clf;
            trisurf(triangles,cx,cy,Z(:,i));
            shading interp;

            if (mode < 4)
                hold on;
                patch('Vertices',[cx cy Z(:,i)],'Faces',rectangles,'FaceVertexCData',Z(:,i),...
                'FaceColor','none');
                hold off;
            end
            
            if (i > 1)
                title(sprintf('nit: %d - error: %e',i-1,err(i)));
            end
            
            if (mode == 2 || mode == 3)
                hold on;
                scatter3(cx,cy,Z(:,i),200/sqrt(size(Z(:,i),1)),'k');
                hold off;
            end

            zlim([zmin zmax]);
            grid on;
            view(vm);
            drawnow;
            pause(0.000001);
        end
        colorbar;
    end
end
