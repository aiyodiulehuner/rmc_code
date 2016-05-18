function sv=NNP_LR_SP(mu,sv,par)
% Computes projection onto mu nuclear norm ball of global matrix X+spZ

global X spZ
[d1,d2]=size(spZ);

while (1)
    [U,S,V]=lansvd('Axz','Atxz',d1,d2,sv,'L');
    diagS = diag(S);
    svp = max(length(find(diagS > 1e-15)),1);
    diagS = diagS(1:svp);
    if svp==par.maxrank
        break
    end
    if svp < sv
        sv = min(svp + 1, par.maxrank);
        break
    else
        sv = min(2*svp, par.maxrank);
    end
    
end
diagS = ProjectOntoL1Ball(diagS, mu,2);
X.U = U(:, 1:svp) * diag(sqrt(diagS(1:svp)));
X.V = V(:, 1:svp) * diag(sqrt(diagS(1:svp)));
