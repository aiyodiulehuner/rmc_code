addpath('NNLS-0/solver');
addpath('NNLS-0/PROPACKmod/');
clear;clc;clear global

%parameters
par.tol     = 1e-3;
par.maxiter = 50;
par.verbose = 1;
f={'spearman_rho', 'kendall_tau', 'NDCG', 'MSE'}%, 'Precision': doesnt make sense without clear thresholds

niter=1:5;

par.maxrank=100;
probiter=0.2:0.2:0.8;
par.nnp = 1;
RMC=1;
SMC=0;
K=10; %NDCG at 10
th=0.5;
if (RMC)
muiter=[5e4,1e4,5000,100,500,100,50,10];
resultRMC=zeros(lenth(niter),length(probiter), length(muiter), 3, length(f));
mu0=1;%sum(svd(Y));
for i=1:length(niter)
    for pi=1:length(probiter)
        p=probiter(pi);
        load(sprintf('../data/neurosynth_counts/folds%d/neurosynth_%d.mat',i,round(p*100)));
        par.maxrank = min([d1,d2,par.maxrank]);
        fprintf('Size: %dX%d, p:%f, train:val:test::%d:%d:%d\n',d1,d2,p,...
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
            [Yrmc,Yrt,iter,res,ii]=rmc_fixed_margin(ii,Jcol,jj,yy,d1,d2,mu,par,Yrmc,yy); 
            t=toc;
            yest=Amap_MatComp(Yrmc,ii,Jcol);
            k1=evalRanking(yy,yest,Jcol,f,K,th);            
            resultRMC(i,pi,m,1,:)=k1;
            
            %validation            
            yest_val=Amap_MatComp(Yrmc,ii_val,Jcol_val);
            k2=evalRanking(yy_val,yest_val,Jcol_val,f,K,th);
            resultRMC(i,pi,m,2,:)=k2;
                        
            %testing
            yest_test=Amap_MatComp(Yrmc,ii_test,Jcol_test);
            k3=evalRanking(yy_test,yest_test,Jcol_test,f,K,th);                
            resultRMC(i,pi,m,3,:)= k3;  
            
            fprintf('RMC  mu:%f. iter:%d, res:%f, ||X||_*:%f, t:%f\n', ...
                mu, iter,res,sum(sum(Yrmc.U.^2)),t);
            save('results/ns_resultRMC.mat','resultRMC')
            save(sprintf('results/ns_yrmc_p%d/ns_yrmc_cv%d_mi%d.mat',round(p*100),i,m),'Yrmc','t')
        end
    end          
end 
end

if SMC
muiter=[5e4,2.5e4,1e4,7500,5000,2500,1000,500,100,50];
resultSMC=zeros(length(niter),length(probiter), length(muiter), 3, length(f));
mu0=1;
niter=1:5;
for i=1:length(niter)
    for pi=1:length(probiter)
        p=probiter(pi);
        load(sprintf('../data/neurosynth_counts/folds%d/neurosynth_%d.mat',i,round(p*100)));
        par.maxrank = min([d1,d2,par.maxrank]);
        fprintf('Size: %dX%d, p:%f, train:val:test::%d:%d:%d\n',d1,d2,p,...
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
            resultSMC(i,pi,m,1,:)=k1;
            
            %validation            
            yest_val=Amap_MatComp(Ysmc,ii_val,Jcol_val);
            k2=evalRanking(yy_val,yest_val,Jcol_val,f,K,th);
            resultSMC(i,pi,m,2,:)=k2;
                        
            %testing
            yest_test=Amap_MatComp(Ysmc,ii_test,Jcol_test);
            k3=evalRanking(yy_test,yest_test,Jcol_test,f,K,th);                
            resultSMC(i,pi,m,3,:)= k3;  
            
            fprintf('SMC  mu:%f. iter:%d, res:%f, ||X||_*:%f, t:%f\n', ...
                mu, iter,res,sum(sum(Ysmc.U.^2)),t);
            save(sprintf('results/ns_resultSMC.mat'),'resultSMC')
            save(sprintf('results/ns_ysmc_p%d/ysmc_cv%d_mi%d.mat',round(p*100),i,m),'Ysmc','t')
        end
    end          
end 
end

