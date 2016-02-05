clear;clc
addpath('pav')
addpath('utils/');
addpath('NNLS-0/solver');
addpath('NNLS-0/PROPACKmod/');
clear global
Y=importdata('../neurosynth_counts/count_matrix.csv')';
[d1,d2]=size(Y);
%parameters
par.tol     = 1e-5;
par.maxiter = 1000;
par.maxrank = 200;%min([d1,d2,500]);
par.verbose = 1;
par.nnp=1;
f={'spearman_rho', 'kendall_tau', 'NDCG'};

niter=1;
muiter=[500,1000,1e4,5e4];
probiter=0.2:0.2:0.8;    
%probiter=[0.4];
resultSMC=zeros(niter,length(probiter), length(muiter), length(f));
resultRMC=zeros(niter,length(probiter), length(muiter), length(f));

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
            %[Ysmc,iter,res]=smc(ii,Jcol,jj,YOmega,d1,d2,mu,par);            
            %k=evalRanking(Y,Ysmc.U*Ysmc.V',f);
            %fprintf('\t SMC mu:%f. iter:%d, res:%f, ||X||_*:%f, ktau:%f, srho:%f, ndcg:%f\n',...
            %    mu, iter,res,sum(sum(Ysmc.U.^2)),k(1),k(2),k(3));
            %resultSMC(i,pi,m,:)= k;      
            [Yrmc,iter,res]=rmc_fixed_margin(ii,Jcol,jj,YOmega,d1,d2,mu,par);
            k=evalRanking(Y,Yrmc.U*Yrmc.V',f);
            fprintf('\t RMC mu:%f. iter:%d, res:%f, ||X||_*:%f, ktau:%f, srho:%f, ndcg:%f\n',...
                mu, iter,res,sum(sum(Yrmc.U.^2)),k(1),k(2),k(3));
            resultRMC(i,pi,m,:)= k;  
        end

    end          
end    
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