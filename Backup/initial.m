clc;clear;clear global;
addpath('utils/');
addpath('NNLS-0/solver');
addpath('NNLS-0/PROPACKmod/');
addpath('pav');
rng(42)

outfile='resultTest.m';

niter=1;
d1=50;d2=50;r=2;
g=@(a)a;
f={'spearman_rho', 'kendall_tau', @(x,y)norm(x-y)^2/length(x-y)};

a=0;b=1;

U=randn(d1,r);
V=randn(d2,r);
theta=U*V'; 
theta=sqrt(d1*d2)*theta/norm(theta,'fro');
Omega=rand(size(theta));

Y=theta;
for j=1:d2
    Y(:,j)=g(theta(:,j));
end

p=1;
[ii,jj]=find(Omega<=p);
YOmega=Y(Omega<=p);
[YOmega,ii,Jcol]=processInput(ii,jj,YOmega);

par.tol     = 1e-2;
par.maxiter = 200;
par.maxrank = min([d1,d2,500]);
par.verbose = 1;

probiter=0.05:0.05:0.95;    

if a
    %muiter=[0.1,0.2,0.5,0.7,1,5,10,20];
    %u0=sum(svd(theta));
    
    muiter=[0.05,0.1,0.5,1,10];
    mu0=1;
    resultSMC=zeros(length(muiter),3);
    for m=1:length(muiter)
        mu=mu0*muiter(m);
        [Ysmc,iter,~]=smc(ii,Jcol,jj,YOmega,d1,d2,mu,par);
        k=evalRanking(theta,Ysmc.U*Ysmc.V',f);k(3)=sqrt(k(3));
        fprintf('Size: %dX%d, rk:%d, p:%f, gi:%s, mu:%f, \n\t iter:%d ktau:%f, srho:%f, rmse:%f, outrank:%d\n',...
            d1,d2,r,p,func2str(g),mu,iter,k(1),k(2),k(3),size(Ysmc.U,2));    
        resultSMC(m,:)=3;
    end
end

if b
    muiter=[0.01];
    mu0=1;
    mu=muiter(1)*mu0;
    resultRMC=zeros(length(muiter),3);
    par.PAV_QP=qpparams_pav_margin(Jcol);
    
    for m=1:length(muiter)
        mu=mu0*muiter(m);
        [Yrmc,Yrt,eps,iter,~,result]=rmc(ii,Jcol,jj,YOmega,d1,d2,mu,par,theta,f);


        k=evalRanking(theta,Yrmc.U*Yrmc.V',f);k(3)=sqrt(k(3));
        fprintf('Size: %dX%d, rk:%d, p:%f, gi:%s, mu:%f, \n\t iter:%d ktau:%f, srho:%f, rmse:%f\n',...
            d1,d2,r,p,func2str(g),mu,iter,k(1),k(2),k(3));    
        resultRMC(m,:)=3;
    end
end