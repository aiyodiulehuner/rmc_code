clear;clc
addpath('utils/');
addpath('NNLS-0/solver');
addpath('NNLS-0/PROPACKmod/');
addpath('pav/')
Y=importdata('../neurosynth_counts/count_matrix.csv')';
[d1,d2]=size(Y);
%parameters
par.tol     = 1e-5;
par.maxiter = 1000;
par.maxrank = 200;%min([d1,d2,500]);
par.verbose = 0;
par.nnp=1;
f={'spearman_rho', 'kendall_tau', 'NDCG'};

niter=1;
muiter=[1e-6,1e-4,1e-3,0.01,0.1,1,10];
probiter=0.2:0.2:1;    
%probiter=[0.4];
resultSMC=zeros(niter,length(probiter), length(muiter), length(f));

Omega=rand(size(Y));

mu0=1;%sum(svd(Y));

p=probiter(1);
[ii,jj]=find(Omega<=p);
        YOmega=Y(Omega<=p);
        [YOmega,ii,Jcol]=processInput(ii,jj,YOmega);
        fprintf('Size: %dX%d, p:%f\n',d1,d2,p);
        mu0=500;
Amap  = @(X) Amap_MatComp(X,ii,Jcol);  
if (length(YOmega)/(d1*d2)>0.6)
    ATmap = @(y) full(sparse(ii,jj,y, d1,d2));
else
    if (exist('mexspconvert')==3); 
        ATmap = @(y) mexspconvert(d1,d2,y,ii,Jcol); 
    else
        ATmap = @(y) sparse(ii,jj,y, d1,d2); 
    end
end

%% Initialize Variables
sv=20; 

global X spZ
rinit=10;
eps=zeros(d2,1);
for j=1:length(Jcol)-1
    ind = Jcol(j)+1:Jcol(j+1);
    eps(j) = max(1e-10,min(diff(YOmega(ind))));
end
%eps=max(eps,1e-10)*ones(d2,1);
Yrt=YOmega;

X.U=zeros(d1,rinit);X.V=zeros(d2,rinit); 
XOmega=Amap(X);
spZ=ATmap((Yrt-XOmega)/2);
Xold=XOmega;
par.continuation=0.5;mu0=mu0/par.continuation;
ch=0; res=0; mu=mu0;
mu=par.continuation*mu;    
sv=NNP_LR_SP(mu,min(sv,par.maxrank));
XOmega=Amap(X);
y=((Yrt+XOmega)/2)';