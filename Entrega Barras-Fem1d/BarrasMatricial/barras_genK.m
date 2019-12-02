function localK = barras_genK(nodes,propiedades)
    
    E = propiedades.E;
    A = propiedades.A;
    
    nodo1 = nodes(1,:);
    nodo2 = nodes(2,:);
    
    alpha = atan2(nodo2(2) - nodo1(2), nodo2(1) - nodo1(1));
    L = norm(nodo2 - nodo1);
    
    c = cos(alpha);
    s = sin(alpha);
    sc = s * c;
    
    k = E*A/L * [c^2 sc; sc s^2];
    localK = [k -k; -k k];
end
