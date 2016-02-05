# rmc_code
addpath('NNLS-0/solver');
addpath('NNLS-0/PROPACKmod/');
addpath('pav');
addpath('utils');

utils/
[yys,iis,Jcol]=processInput(ii,jj,yy)
returns csc format: iis are indices, Jcol is indptr (starts at 0
ii(Jcol(j)+1:Jcol(j+1)) is ordered such that yy(Jcol(j)+1:Jcol(j+1)) is sorted min->max

evalRanking:
implementation of Spearman Rho, Kendall Tau, NDCG@k

Atxz, Axz 
A=X.U*X.V'+spZ X,spZ are global variables
returns A'z, Az


NNLS-0/solver/
x=Amap_MatComp(X,ii,Jcol) -> X->POmega(X)


