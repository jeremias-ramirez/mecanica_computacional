function [reaction] = fem2d_pstr_reaction(K,F,U,Fixnodes)
    if isempty(Fixnodes)
        disp('************');
        disp('El problema ha sido mal condicionado. Se espera al menos una condición Dirichlet');
        disp('************');
        reaction = [];
    else
        reaction = sparse(size(U,1),1);
        fix = 2*Fixnodes(:,1) + Fixnodes(:,2) - 2;
        reaction(fix) = K(fix,:) * U - F(fix);
    end
end