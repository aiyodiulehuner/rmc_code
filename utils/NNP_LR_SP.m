function sv=NNP_LR_SP(mu,sv)
% Computes projection onto mu nuclear norm ball of global matrix X+spZ

global X spZ
[d1,d2]=size(spZ);

[X.U,S,X.V]=lansvd('Axz','Atxz',d1,d2,sv,'L');

diagS = diag(S);
svp = max(length(find(diagS>1e-15)),1);
diagS = diagS(1:svp);

diagS = ProjectOntoL1Ball(diagS, mu,2);
svp=max(length(find(diagS>1e-15)),1);
if svp < sv %|| iter < 10
    sv = min(svp + 1, d2);
    %sv=sv;
else
    sv = min(2*svp, d2);
end

X.U = X.U(:, 1:svp) * diag(sqrt(diagS(1:svp)));
X.V = X.V(:, 1:svp) * diag(sqrt(diagS(1:svp)));