function [sv,svp]=SVT_LR_SP(th,sv,options)
% Computes SVT of global matrix X+spZ by threhold th

global X spZ
[d1,d2]=size(spZ);

[X.U,S,X.V]=lansvd('Axz','Atxz',d1,d2,sv,'L',options);
% predict sv, ref: Toh et al. 2009,
% In some iterations this is not exactly mu SVT: if S(svp+1)>=mu
diagS = diag(S);
diagS = diagS(1:sv);svn = length(find(diagS > th));
diagS = diagS(1:max(svn,1));
svp = max(svn,1);
%if (svp>1)
%    ratio = diagS(1:end-1)./diagS(2:end);
%    [max_ratio, max_idx] = max(ratio);
%    if max_ratio > 5
%        svp = min(svn, max_idx);
%    end
%end
if svp < sv %|| iter < 10
    sv = min(svp + 1, d2);
else
    sv = min(svp + 10, d2);
end
if (svp>1)
    sqrtds = sqrt(max(0,diagS(1:svp) - th));
else
    sqrtds = sqrt(diagS(1:svp));
end
X.U = X.U(:, 1:svp) * diag(sqrtds);
X.V = X.V(:, 1:svp) * diag(sqrtds);