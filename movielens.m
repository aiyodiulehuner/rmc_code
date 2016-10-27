clear;clc
addpath('pav')
addpath('utils/');
addpath('NNLS-0/solver');
addpath('NNLS-0/PROPACKmod/');
clear global

%parameters
par.tol     = 1e-3;
par.maxiter = 50;
par.verbose = 1;
f={'spearman_rho', 'kendall_tau', 'NDCG', 'MSE', 'Precision'};


par.maxrank=1000;
cviter=[1,2,3,4,5];
par.nnp = 1;
RMC=0;
SMC=1;
K=5;
th=3;
if (RMC)
    muiter=[1,5,10,20,50,100,250,500,1000];    
    resultRMC=zeros(cviter, length(muiter), 3, length(f));
    mu0=1;%sum(svd(Y));

    for ci=1:length(cviter)
        cv=cviter(ci);
        load(sprintf('../data/ml-100k/folds/ml_%d.mat',ci));
        par.maxrank = min([d1,d2,par.maxrank]);
        fprintf('Size: %dX%d, p:%f, train:val:test::%d:%d:%d\n',d1,d2,cv,...
            length(yy),length(yy_val),length(yy_test));
        %d1,d2,yy,ii,Jcol,yy_val,ii_val,Jcol_val,yy_test,ii_test,Jcol_test
        ii_train=ii;
        Yrmc.U=zeros(d1,10);Yrmc.V=zeros(d2,10);
        Yrt=yy;
        for m=1:length(muiter)           
            mu=mu0*muiter(m);
            fprintf('mu=%f, nnp:%d\n',mu,par.nnp)
            % training
            tic;
            [Yrmc,Yrt,iter,res,ii]=rmc_fixed_margin_AM(ii,Jcol,jj,yy,d1,d2,mu,par,Yrmc,Yrt); 
            t=toc;
            yest=Amap_MatComp(Yrmc,ii,Jcol);            
            k1=evalRanking(yy,yest,Jcol,f,K,th);            
            resultRMC(ci,m,1,:)=k1;
            
            %validation            
            yest_val=Amap_MatComp(Yrmc,ii_val,Jcol_val);
            k2=evalRanking(yy_val,yest_val,Jcol_val,f,K,th);
            resultRMC(ci,m,2,:)=k2;
                        
            %testing
            yest_test=Amap_MatComp(Yrmc,ii_test,Jcol_test);
            k3=evalRanking(yy_test,yest_test,Jcol_test,f,K,th);                
            resultRMC(ci,m,3,:)= k3;  
            
            fprintf('RMC  mu:%f. iter:%d, res:%f, ||X||_*:%f, t:%f\n', ...
                mu, iter,res,sum(sum(Yrmc.U.^2)),t);
            save('results/ml_resultRMC.mat','resultRMC')
            save(sprintf('results/ml_yrmc_cv%d_mi%d.mat',ci,m),'Yrmc','t')
        end
    end          
end 

if SMC
    muiter=[5e4,2.8e4,1e4,5000,1000,500,100,50];
    resultSMC=zeros(length(cviter), length(muiter), 3, length(f));
    mu0=1;%sum(svd(Y));

    for ci=1:length(cviter)
        cv=cviter(ci);
        load(sprintf('../data/ml-100k/folds/ml_%d.mat',ci));
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
            resultSMC(ci,m,1,:)=k1;
            
            %validation            
            yest_val=Amap_MatComp(Ysmc,ii_val,Jcol_val);
            k2=evalRanking(yy_val,yest_val,Jcol_val,f,K,th);
            resultSMC(ci,m,2,:)=k2;
                        
            %testing
            yest_test=Amap_MatComp(Ysmc,ii_test,Jcol_test);
            k3=evalRanking(yy_test,yest_test,Jcol_test,f,K,th);                
            resultSMC(ci,m,3,:)= k3;  
            
            fprintf('SMC  mu:%f. iter:%d, res:%f, ||X||_*:%f, t:%f\n', ...
                mu, iter,res,sum(sum(Ysmc.U.^2)),t);
            save('results/ml_resultSMC.mat','resultSMC')
            save(sprintf('results/ml_ysmc_cv%d_mi%d.mat',ci,m),'Ysmc','t')
        end
    end          
end 
