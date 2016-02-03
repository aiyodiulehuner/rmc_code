clear;clear global;
addpath('utils/');
addpath('NNLS-0/solver');
addpath('NNLS-0/PROPACKmod/');
addpath('pav');

rmc = 1;
niter = 1;
d1 = 50; d2 = 50; r = 2;
%U = randn(d1,r);
%V = randn(d2,r);
%theta = U*V'; 
%Y = sqrt(d1*d2)*theta/norm(theta,'fro');
load('theta'); load('Y'); load('Yrmc');

f = {'spearman_rho', 'kendall_tau', @(x,y)norm(x-y)^2/length(x-y)};
outfile = 'resultTest.m';
piter = 0.9;
result = zeros(niter,length(piter),length(f)+1);



par.tol     = 1e-6;
par.maxiter = 500;
par.maxrank = min([d1,d2,500]);
par.verbose = 0;
par.continuation = 0.5;
par.mu0 = 10;%3.306146
%par.mutarget=3.3;

debug = {theta,f};

Yold = zeros(size(theta));

for n=1:niter
    rng('shuffle')
    Omega=rand(size(theta));    
    for pi=1:length(piter)
        p=piter(pi);
        [ii,jj]=find(Omega<=p);
        YOmega=Y(Omega<=p);
        [YOmega,ii,Jcol]=processInput(ii,jj,YOmega);
        YrOmega=Amap_MatComp(Yrmc,ii,Jcol); 
        eps=zeros(d2,1);
        for j=1:length(Jcol)-1
            ind = Jcol(j)+1:Jcol(j+1);
            if length(ind)<2
                continue
            end
            eps(j) = min(diff(YOmega(ind)));        
        end        
        assert(all(eps>0))
        %eps=ones(size(eps))*(1.0/d2);
        [X,spZ]=RMC_exact_fixed_margin(ii,jj,Jcol,YOmega,eps,d1,d2,par,debug);
    end
    %% Results
    if ~isempty(debug)
        Yrmc=X.U*X.V'+full(spZ);
        norm(Yold-Yrmc)
        Yold=Yrmc;
    end
    k=evalRanking(theta,Yrmc,f);k(3)=sqrt(k(3));
    k1=evalRanking(theta,X.U*X.V',f);k(3)=sqrt(k(3));
    s=svd(Yrmc);
    fprintf('Size: %dX%d, rk:%d, p:%f, mu:(%f,%f), normspZ:%f \n\t ktau:(%f,%f), srho:(%f,%f), rmse:(%f,%f)\n',...
        d1,d2,r,p,sum(svd(Yrmc)),sum(svd(X.U*X.V')),norm(spZ,'fro'),k(1),k1(1),k(2),k1(2),k(3),k1(3));    
    result(n,pi,1:3)=k;
    result(n,pi,4)=find(cumsum(s.^2)/sum(s.^2)>(1-1e-5),1); 
end
