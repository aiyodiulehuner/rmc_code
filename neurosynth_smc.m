clear;clc
addpath('pav')
addpath('utils/');
addpath('NNLS-0/solver');
addpath('NNLS-0/PROPACKmod/');
clear global

%parameters
par.tol     = 1e-3;
par.maxiter = 1000;

par.verbose = 1;
par.nnp = 1;
f={'spearman_rho', 'kendall_tau', 'NDCG'};

niter=1;
muiter=[250,500,1000,1e4,5e4];
probiter=0.8;%:0.2:0.8;    
%probiter=[0.4];
resultSMC=zeros(niter,length(probiter), length(muiter), 3, length(f));
resultRMC=zeros(niter,length(probiter), length(muiter), 3, length(f));


mu0=1;%sum(svd(Y));


for i=1:length(niter)
    for pi=1:length(probiter)
        p=probiter(pi);
        load(sprintf('../neurosynth_counts/folds/neurosynth_%d.mat',round(p*100)));
        par.maxrank = min(d1,d2);
        fprintf('Size: %dX%d, p:%f, train:val:test::%d:%d:%d\n',d1,d2,p,...
            length(yy),length(yy_val),length(yy_test));
        %d1,d2,yy,ii,Jcol,yy_val,ii_val,Jcol_val,yy_test,ii_test,Jcol_test
        ii_train=ii;         
        for m=1:length(muiter)           
            mu=mu0*muiter(m);
            fprintf('mu=%f\n',mu)
            % training
            [Yrmc,Yrt,iter,res]=rmc_fixed_margin(ii,Jcol,jj,yy,d1,d2,mu,par);            
            yest=Amap_MatComp(Yrmc,ii_train,Jcol);
            
            k1=evalRanking(yy,yest,Jcol,f);            
            resultRMC(i,pi,m,1,:)=k1;
            
            %validation            
            yest_val=Amap_MatComp(Yrmc,ii_val,Jcol_val);
            k2=evalRanking(yy_val,yest_val,Jcol_val,f);
            resultRMC(i,pi,m,2,:)=k2;
                        
            %testing
            yest_test=Amap_MatComp(Yrmc,ii_test,Jcol_test);
            k3=evalRanking(yy_test,yest_test,Jcol_test,f);                
            resultRMC(i,pi,m,3,:)= k3;  
            
            fprintf('RMC  mu:%f. iter:%d, res:%f, ||X||_*:%f', ...
                mu, iter,res,sum(sum(Yrmc.U.^2)));
            save('resultRMC.mat','resultRMC')
            save(sprintf('yrmc_%d_%d.mat',round(p*100),m),'Yrmc')
        end

    end          
end