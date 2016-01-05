clc;clear;clear global;
addpath('utils/');
addpath('NNLS-0/solver');
addpath('NNLS-0/PROPACKmod/');
addpath('pav');
rng(42)

outfile='resultTest.m';

niter=1;
d1=100;d2=100;r=5;
g=@(a)a;
f={'spearman_rho', 'kendall_tau', @(x,y)norm(x-y)^2/length(x-y)};

U=randn(d1,r);
V=randn(d2,r);
theta=U*V'; 
theta=sqrt(d1*d2)*theta/norm(theta,'fro');
Omega=rand(size(theta));

Y=theta;
for j=1:d2
    Y(:,j)=g(theta(:,j));
end

par.tol     = 1e-3;
par.maxiter = 1000;
par.maxrank = min([d1,d2,500]);
par.verbose = 1;

muiter=[0.05];
mu0=1;
resultRMC=zeros(length(muiter),3);
par.PAV_QP=qpparams(Jcol);
for m=1:length(muiter)
    mu=mu0*muiter(m);
    [Ysmc,iter,~,result]=rmc(ii,Jcol,jj,YOmega,d1,d2,mu,par,theta,f);


    k=evalRanking(theta,Ysmc.U*Ysmc.V',f);k(3)=sqrt(k(3));
    fprintf('Size: %dX%d, rk:%d, p:%f, gi:%s, mu:%f, \n\t iter:%d ktau:%f, srho:%f, rmse:%f\n',...
        d1,d2,r,p,func2str(g),mu,iter,k(1),k(2),k(3));    
    resultRMC(m,:)=3;
end
