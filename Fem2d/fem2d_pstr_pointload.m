function [F] = fem2d_pstr_pointload(F,pointload,xnode,icone)
    for n = 1 : size(pointload,1)
        e = pointload(n,1);
        Fx = pointload(n,2);
        Fy = pointload(n,3);
        xp = pointload(n,4);
        yp = pointload(n,5);
        if icone(e,4) == -1 % triangular element
            ele = icone(e,1:3);
        else
            ele = icone(e,:);
        end
        nodes = xnode(ele,:);
        N = fem2d_pstr_blerp(nodes,xp,yp);
        indx = [];
        f = [];
        for i=1:length(ele)
            indx = [indx ele(i)*2-1 ele(i)*2];
            f = [f N(i)*Fx N(i)*Fy];
        end
        F(indx) = F(indx) + f';
    end
end