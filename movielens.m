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
f={'spearman_rho', 'kendall_tau', 'NDCG', 'MSE', 'Precision'};
niter=1;

par.maxrank=1000;
cviter=1;
par.nnp = 1;
RMC=0;
SMC=1;
K=5;
th=3;
if (RMC)
muiter=[250,500,1000,1e4,5e4];    
resultRMC=zeros(niter,length(cviter), length(muiter), 3, length(f));
mu0=1;%sum(svd(Y));
for i=1:length(niter)
    for ci=1:length(cviter)
        cv=cviter(ci);
        load(sprintf('../ml-100k/folds/ml_a.mat'));
        par.maxrank = min([d1,d2,par.maxrank]);
        fprintf('Size: %dX%d, p:%f, train:val:test::%d:%d:%d\n',d1,d2,cv,...
            length(yy),length(yy_val),length(yy_test));
        %d1,d2,yy,ii,Jcol,yy_val,ii_val,Jcol_val,yy_test,ii_test,Jcol_test
        ii_train=ii;         
        for m=1:length(muiter)           
            mu=mu0*muiter(m);
            fprintf('mu=%f, nnp:%d\n',mu,par.nnp)
            % training
            tic;
            [Yrmc,Yrt,iter,res]=rmc_fixed_margin(ii,Jcol,jj,yy,d1,d2,mu,par);            
            t=toc;
            yest=Amap_MatComp(Yrmc,ii_train,Jcol);            
            k1=evalRanking(yy,yest,Jcol,f,K,th);            
            resultRMC(i,ci,m,1,:)=k1;
            
            %validation            
            yest_val=Amap_MatComp(Yrmc,ii_val,Jcol_val);
            k2=evalRanking(yy_val,yest_val,Jcol_val,f,K,th);
            resultRMC(i,ci,m,2,:)=k2;
                        
            %testing
            yest_test=Amap_MatComp(Yrmc,ii_test,Jcol_test);
            k3=evalRanking(yy_test,yest_test,Jcol_test,f,K,th);                
            resultRMC(i,ci,m,3,:)= k3;  
            
            fprintf('RMC  mu:%f. iter:%d, res:%f, ||X||_*:%f, t:%f\n', ...
                mu, iter,res,sum(sum(Yrmc.U.^2)),t);
            save('ml_resultRMC.mat','resultRMC')
            save(sprintf('ml_yrmc_a_%d.mat',m),'Yrmc','t')
        end
    end          
end 
end

if SMC
muiter=[5e4,2.8e4,1e4,5000,1000,500,100,50];
resultSMC=zeros(niter,length(cviter), length(muiter), 3, length(f));
mu0=1;%sum(svd(Y));

for i=1:length(niter)
    for ci=1:length(cviter)
        cv=cviter(ci);
        load(sprintf('../ml-100k/folds/ml_a.mat'));
        par.maxrank = min([d1,d2,par.maxrank]);
        fprintf('Size: %dX%d, p:%f, train:val:test::%d:%d:%d\n',d1,d2,cv,...
            length(yy),length(yy_val),length(yy_test));
        %d1,d2,yy,ii,Jcol,yy_val,ii_val,Jcol_val,yy_test,ii_test,Jcol_test
        ii_train=ii;         
        Ysmc.U=zeros(d1,10);Ysmc.V=zeros(d2,10);
        for m=1:length(muiter)           
            mu=mu0*muiter(m);
            fprintf('mu=%f, nnp:%d\n',mu,par.nnp)
            % training
            tic;
            [Ysmc,iter,res]=smc(ii,Jcol,jj,yy,d1,d2,mu,par,Ysmc);            
            t=toc;
            yest=Amap_MatComp(Ysmc,ii_train,Jcol);            
            k1=evalRanking(yy,yest,Jcol,f,K,th);            
            resultSMC(i,ci,m,1,:)=k1;
            
            %validation            
            yest_val=Amap_MatComp(Ysmc,ii_val,Jcol_val);
            k2=evalRanking(yy_val,yest_val,Jcol_val,f,K,th);
            resultSMC(i,ci,m,2,:)=k2;
                        
            %testing
            yest_test=Amap_MatComp(Ysmc,ii_test,Jcol_test);
            k3=evalRanking(yy_test,yest_test,Jcol_test,f,K,th);                
            resultSMC(i,ci,m,3,:)= k3;  
            
            fprintf('SMC  mu:%f. iter:%d, res:%f, ||X||_*:%f, t:%f\n', ...
                mu, iter,res,sum(sum(Ysmc.U.^2)),t);
            save('ml_resultSMC.mat','resultSMC')
            save(sprintf('ml_ysmc_a_%d.mat',m),'Ysmc','t')
        end
    end          
end 
end