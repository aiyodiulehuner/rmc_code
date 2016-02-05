clear;clc
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

for i=1:length(niter)
    for pi=1:length(probiter)
        p=probiter(pi);

        [ii,jj]=find(Omega<=p);
        YOmega=Y(Omega<=p);
        [YOmega,ii,Jcol]=processInput(ii,jj,YOmega);
        fprintf('Size: %dX%d, p:%f\n',d1,d2,p);

        for m=1:length(muiter)
            mu=mu0*muiter(m);
            [Ysmc,iter,res]=smc(ii,Jcol,jj,YOmega,d1,d2,mu,par);
            k=evalRanking(Y,Ysmc.U*Ysmc.V',f);
            fprintf('\t mu:%f. iter:%d, res:%f, ||X||_*:%f, ktau:%f, srho:%f, ndcg:%f\n',...
                mu, iter,res,sum(sum(Ysmc.U.^2)),k(1),k(2),k(3));
            resultSMC(i,pi,m,:)= k;                
        end

    end          
end    
