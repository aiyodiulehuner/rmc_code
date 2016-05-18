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
probiter=[0.8];%:0.2:0.8;
par.nnp = 1;
RMC=1;
SMC=0;

if (RMC)
muiter=[5e4,1e4,5000,1000,500,250,100,50];
resultRMC=zeros(niter,length(probiter), length(muiter), 3, length(f));
mu0=1;%sum(svd(Y));
for i=1:length(niter)
    for pi=1:length(probiter)
        p=probiter(pi);
        load(sprintf('../neurosynth_counts/folds/neurosynth_%d.mat',round(p*100)));
        par.maxrank = min([d1,d2,par.maxrank]);
        fprintf('Size: %dX%d, p:%f, train:val:test::%d:%d:%d\n',d1,d2,p,...
            length(yy),length(yy_val),length(yy_test));
        %d1,d2,yy,ii,Jcol,yy_val,ii_val,Jcol_val,yy_test,ii_test,Jcol_test
        ii_train=ii;
        Yrmc.U=zeros(d1,10);Yrmc.V=zeros(d2,10);
        for m=1:length(muiter)           
            mu=mu0*muiter(m);
            fprintf('mu=%f, nnp:%d\n',mu,par.nnp)
            % training
            tic;
            [Yrmc,Yrt,iter,res]=rmc_fixed_margin(ii,Jcol,jj,yy,d1,d2,mu,par,Yrmc); 
            t=toc;
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
            
            fprintf('RMC  mu:%f. iter:%d, res:%f, ||X||_*:%f, t:%f\n', ...
                mu, iter,res,sum(sum(Yrmc.U.^2)),t);
            save('resultRMC.mat','resultRMC')
            save(sprintf('yrmc_%d_%d.mat',round(p*100),m),'Yrmc','t')
        end
    end          
end 
end

if SMC
muiter=[5e4,1e4,5000,1000,500,250,100,50];  
resultSMC=zeros(niter,length(probiter), length(muiter), 3, length(f));
mu0=1;%sum(svd(Y));
for i=1:length(niter)
    for pi=1:length(probiter)
        p=probiter(pi);
        load(sprintf('../neurosynth_counts/folds/neurosynth_%d.mat',round(p*100)));
        par.maxrank = min([d1,d2,par.maxrank]);
        fprintf('Size: %dX%d, p:%f, train:val:test::%d:%d:%d\n',d1,d2,p,...
            length(yy),length(yy_val),length(yy_test));
        %d1,d2,yy,ii,Jcol,yy_val,ii_val,Jcol_val,yy_test,ii_test,Jcol_test
        ii_train=ii;         
        Ysmc.U=zeros(d1,10);Ysmc.V=zeros(d2,10)
        for m=1:length(muiter)           
            mu=mu0*muiter(m);
            fprintf('mu=%f, nnp:%d\n',mu,par.nnp)
            % training
            tic;
            [Ysmc,iter,res]=smc(ii,Jcol,jj,yy,d1,d2,mu,par,Ysmc);            
            t=toc;
            yest=Amap_MatComp(Ysmc,ii_train,Jcol);            
            k1=evalRanking(yy,yest,Jcol,f);            
            resultSMC(i,pi,m,1,:)=k1;
            
            %validation            
            yest_val=Amap_MatComp(Ysmc,ii_val,Jcol_val);
            k2=evalRanking(yy_val,yest_val,Jcol_val,f);
            resultSMC(i,pi,m,2,:)=k2;
                        
            %testing
            yest_test=Amap_MatComp(Ysmc,ii_test,Jcol_test);
            k3=evalRanking(yy_test,yest_test,Jcol_test,f);                
            resultSMC(i,pi,m,3,:)= k3;  
            
            fprintf('SMC  mu:%f. iter:%d, res:%f, ||X||_*:%f, t:%f\n', ...
                mu, iter,res,sum(sum(Ysmc.U.^2)),t);
            save('resultSMC.mat','resultSMC')
            save(sprintf('ysmc_%d_%d.mat',round(p*100),m),'Ysmc','t')
        end
    end          
end 
end






plt=0;
if (plt)
    load('resultSMC')
    resultbestSMC=zeros(niter, length(probiter),length(f));
    for i=1:niter
        for pi=1:length(probiter)
            temp=squeeze(resultSMC(i,pi,:,:));
            [k,m]=max(temp);
            %k(3)=k(3)^2;
            resultbestSMC(i,pi,:)= k;
        end
    end
    resultbestRMC=zeros(niter, length(probiter),length(f));
    for i=1:niter
        for pi=1:length(probiter)
            temp=squeeze(resultRMC(i,pi,:,:));
            [k,m]=max(temp);
            %k(3)=k(3)^2;
            resultbestRMC(i,pi,:)= k;
        end
    end
    f={'Spearman \rho', 'Kendall \tau', 'NDCG\@50'};
    figure()
    for k=1:length(f)
        subplot(1,length(f),k);
        plot(probiter,squeeze(resultbestSMC(i,:,k)),'-o','LineWidth',2); hold on;
        plot(probiter,squeeze(resultbestRMC(i,:,k)),'-o','LineWidth',2); hold on;
        xlabel('$\frac{|\Omega|}{2d_1d_2}$','Interpreter','latex','FontSize',16)
        ylabel(f(k),'FontSize',16)           
    end
    legend({'SMC','RMC'})
end
