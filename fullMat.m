clc;clear;clear global;
addpath('utils/');
addpath('NNLS-0/solver');
addpath('NNLS-0/PROPACKmod/');
addpath('pav');
rng(42)

outfile='resultTest.m';

niter=1;
d1=50;d2=50;r=2;
f={'spearman_rho', 'kendall_tau', @(x,y)norm(x-y)^2/length(x-y)};

U=randn(d1,r);
V=randn(d2,r);
theta=U*V'; 
Y=sqrt(d1*d2)*theta/norm(theta,'fro');
Omega=rand(size(theta));

par.tol     = 1e-3;
par.maxiter = 1000;
par.maxrank = min([d1,d2,500]);
par.verbose = 1;

p=1;
[ii,jj]=find(Omega<=p);
YOmega=Y(Omega<=p);
[YOmega,ii,Jcol]=processInput(ii,jj,YOmega);



%% RMC 
Amap  = @(X) Amap_MatComp(X,ii,Jcol);  
if (length(YOmega)/(d1*d2)>0.6)
    ATmap = @(y) full(sparse(ii,jj,y, d1,d2));
else
    if (exist('mexspconvert')==3); 
        ATmap = @(y) mexspconvert(d1,d2,y,ii,Jcol); 
    else
        ATmap = @(y) sparse(ii,jj,y, d1,d2); 
    end
end


global X spZ
sv=20;
t=1;%stepsize

eps=zeros(d2,1);
for j=1:length(Jcol)-1
    ind = Jcol(j)+1:Jcol(j+1);
    eps(j) = min(diff(YOmega(ind)));
end
Yrt=YOmega;

%% Assign Variables
rinit=10;X.U=zeros(d1,rinit);X.V=zeros(d2,rinit);
XOmega=Amap(X);
spZ=ATmap(t*(Yrt-XOmega)); 

mu=1; 
 %% Tracking
XOld=XOmega;YOld=Yrt; % for exit conditions
ff=[0.5*norm(XOmega-Yrt)^2,mu*sum(sum(X.U.^2))];
disp(ff)
result=[];


%% Iterations
for iter=1:par.maxiter
    %% UPDATE                  
    %Xupdate
    spZ=ATmap(t*(Yrt-XOmega));     
    sv=SVT_LR_SP(t*mu,min(sv,par.maxrank));%SVT of X+spZ  
    XOmega=Amap(X);
    %Yupdate      
    Yrt=c_colMR_fixed_margin(XOmega',eps',Jcol'); Yrt=Yrt';


    %% EXIT CONDITIONS
    ch=norm(XOld-XOmega)/sqrt(length(XOmega));         
    ff=[0.5*norm(XOmega-Yrt)^2,mu*sum(sum(X.U.^2))];
    if par.verbose
        fprintf('iter:%d,sv:%d,eps:%2g, ch1:%f, ch2:%0.2g,fval:%0.5g:%0.5g:%0.5g\n',...
            iter,sv,min(eps),ch,norm(Yrt-YOld),sum(ff),ff(1),ff(2))
        k=evalRanking(theta,X.U*X.V',f);k(3)=sqrt(k(3));%cmt
        result=[result;[k',min(eps)]];%cmt
    end     
    XOld=XOmega;YOld=Yrt;
    if ((ch<par.tol) || (iter==par.maxiter))          
        break
    end                

end

Yrmc=X;
k=evalRanking(theta,Yrmc.U*Yrmc.V',f);k(3)=sqrt(k(3));
fprintf('Size: %dX%d, rk:%d, p:%f, mu:%f, \n\t iter:%d ktau:%f, srho:%f, rmse:%f\n',...
        d1,d2,r,p,mu,iter,k(1),k(2),k(3));    
