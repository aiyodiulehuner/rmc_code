%% Data pre-processing
addpath('NNLS-0/solver');
addpath('NNLS-0/PROPACKmod/');
addpath('NNLS-0/PROPACKmod/');
close all;
clc;clear;clear global;

d=[500,500]; rk=5;
g='argsort';
%muiter=[0.01,0.05,0.1,0.2,0.5,1];

muiter=[0.01,0.05,0.1,0.2,0.5,0.7,1,5,10,20];
niter=1;
f={@(x,y)corr(x,y,'type','Kendall'),@(x,y)corr(x,y,'type','spearman')};
probiter=(0.1:0.1)*(sum(d)*rk*log(sum(d)))/(d(1)*d(2));
evalMetrics=zeros(niter,length(probiter),length(muiter),length(f));

for it=1:niter
    data=simulatedData(d,rk);
    %mu=data.mu;
    theta=data.U*data.V';
   
    if (strcmp(g,'I')) 
        mu0=data.mu;
    elseif (strcmp(g,'exp')) 
        mu0=sum(svd(exp(theta)));        
    elseif (strcmp(g,'log')) 
        m=min(min(theta));
        mu0=1*sum(svd(log(theta-m+1)));
    elseif (strcmp(g,'argsort'))
        Ytemp=theta;
        for j=1:d(2)
            Yj=Ytemp(:,j);
            [y,IX]=sort(Yj,'ascend');
            i=1:length(IX);
            Ytemp(:,j)=i(IX);
        end
        mu0=0.1*sum(svd(Ytemp));
        clear Ytemp
    end
    for p=1:length(probiter)
        prob=probiter(p);
        fprintf('Size:%dX%d, rk:%d, p:%f, g:%s, mu0:%2g\n', d(1),d(2),rk,prob,g,mu0);

        [ii,jj]=find(data.Omega<=prob);
        YOmega=theta(data.Omega<=prob);
        n=length(YOmega);
        assert(length(ii)==n);assert(issorted(jj));
        Jcol=compJcol(jj); %csc sparse format ii is indices, jcol is indptr
        for j=1:d(2)
            Yj=YOmega(Jcol(j)+1:Jcol(j+1));
            [y,IX]=sort(Yj,'ascend');
            Ij=ii(Jcol(j)+1:Jcol(j+1));
            ii(Jcol(j)+1:Jcol(j+1))=Ij(IX);
            if (strcmp(g,'argsort'))
                YOmega(Jcol(j)+1:Jcol(j+1))=(1:length(IX))-mean(1:length(IX));
            elseif (strcmp(g,'I')) 
                YOmega(Jcol(j)+1:Jcol(j+1))=y;
            elseif (strcmp(g,'exp')) 
                YOmega(Jcol(j)+1:Jcol(j+1))=exp(y);
            elseif (strcmp(g,'log')) 
                m=min(min(theta));
                YOmega(Jcol(j)+1:Jcol(j+1))=log(y-m+1);
            end
        end
        %-----------
        % SMC
        %-----------
        
        par.tol     = 1e-5;
        par.maxiter = 1000;
        par.maxrank = 500;
        par.verbose = 0;
        for m=1:length(muiter)
            mu=mu0*muiter(m);
            fprintf('mu=%2g\n', mu)            
            [result.Yest,iter,res]=smc(ii,Jcol,jj,YOmega,d(1),d(2),mu,par);
            result.k=evalRanking(theta,result.Yest.U*result.Yest.V',f);
            result.p=prob;
            disp(['SMC exit stats: niter=' num2str(iter) ', residual=' num2str(res) ', rk=' num2str(size(result.Yest.U,2))]);
            disp(['Kendal Tau and Spearman Rho: ' num2str(result.k(1)) ' ' num2str(result.k(2))]);
            evalMetrics(it,p,m,:)=result.k;
            %save(sprintf('results_smc_%d_%d_%d.mat', it,p,m),'result');
        end
    end
end