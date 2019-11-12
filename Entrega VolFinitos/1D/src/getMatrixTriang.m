function k = getMatrixTriang(N,a,b,c)
	k = diag(b*ones(N,1))+diag(c*ones(N-1,1),1)+diag(a*ones(N-1,1),-1);
end

