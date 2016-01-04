clc;clear;clear global;
addpath('utils/');
addpath('NNLS-0/solver');
addpath('NNLS-0/PROPACKmod/');
addpath('pav');



outfile='resultTest.m';

niter=1;
d1=100;d2=100;r=5;
g=@(a)a;

U=randn(d1,r);
V=randn(d2,r);
theta=U*V'; 
theta=sqrt(d1*d2)*theta/norm(theta,'fro');

mu0=sum(svd(theta));%Check
Omega=rand(size(theta));


Y=theta;
for j=1:d2
    Y(:,j)=g(theta(:,j));
end

p=0.4;
[ii,jj]=find(Omega<=p);
YOmega=Y(Omega<=p);
[YOmega,ii,Jcol]=processInput(ii,jj,YOmega);

par.tol     = 1e-5;
par.maxiter = 1000;
par.maxrank = min([d1,d2,500]);
par.verbose = 1;

muiter=[0.01,0.05,0.1,0.2,0.5,0.7,1,5,10,20];
mu0=sum(svd(theta));

r=zeros(length(muiter),3);
mu=mu0;
for m=1:length(muiter)
    mu=mu0*muiter(m);
    [Ysmc,iter,~]=smc(ii,Jcol,jj,YOmega,d1,d2,mu,par);
    k=evalRanking(theta,Ysmc.U*Ysmc.V',f);k(3)=sqrt(k(3));
    fprintf('Size: %dX%d, rk:%d, p:%f, gi:%d, mu:%f, \n\t iter:%d ktau:%f, srho:%f, rmse:%f\n',...
        d1,d2,r,p,gi,mu,iter,k(1),k(2),k(3));    
    r(m,:)=3;
end