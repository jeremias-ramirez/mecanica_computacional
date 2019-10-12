function [dt] = fvm2d_explicit_delta_t(cells,k,rho,cp)
    M = size(cells,2);
    alpha = mean(k)/(rho*cp);
    nd = 2;                     % porque estamos en 2D
    
    for i = 1 : M
        dx = cells(i).dx;
        dy = cells(i).dy;
        delta = mean([dx,dy]);
    end

    dt = 0.5*mean(delta)^2/(alpha*nd);
%     dt = mean(delta)/(alpha*nd);
%     dt = 0.5*abs(1/log10(mean(delta)^2))/(alpha*nd);
end

