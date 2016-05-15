function sv=SVT_LR_SP(th,sv,par)
% Computes SVT of global matrix X+spZ by threhold th

global X spZ
[d1,d2]=size(spZ);

while (1)
    [U,S,V]=lansvd('Axz','Atxz',d1,d2,sv,'L');
    diagS = diag(S);
    diagS = diagS(1:sv);svn = length(find(diagS > th));
    diagS = diagS(1:max(svn,1));
    svp = max(svn,1);
    if svp==par.maxrank
        break
    end
    if svp < sv
        sv = min(svp + 1, par.maxrank);
        break
    else
        sv = min(svp+50, par.maxrank);
    end
    
end
sqrtds = sqrt(max(0,diagS(1:svp) - th));
X.U = U(:, 1:svp) * diag(sqrtds);
X.V = V(:, 1:svp) * diag(sqrtds);

% predict sv, ref: Toh et al. 2009,
% In some iterations this is not exactly mu SVT: if S(svp+1)>=mu

%if (svp>1)
%    ratio = diagS(1:end-1)./diagS(2:end);
%    [max_ratio, max_idx] = max(ratio);
%    if max_ratio > 5
%        svp = min(svn, max_idx);
%    end
%end