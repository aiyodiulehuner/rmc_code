clc;clear;clear global;
addpath('utils/');
addpath('NNLS-0/solver');
addpath('NNLS-0/PROPACKmod/');
addpath('pav');
rng('shuffle')

outfile='resultTest.m';
rmc=1;
niter=1;
d1=50;d2=50;r=2;
f={'spearman_rho', 'kendall_tau', @(x,y)norm(x-y)^2/length(x-y)};

U=randn(d1,r);
V=randn(d2,r);
theta=U*V'; 
Y=sqrt(d1*d2)*theta/norm(theta,'fro');
Omega=rand(size(theta));

par.tol     = 1e-5;
par.maxiter = 1000;
par.maxrank = min([d1,d2,500]);
par.verbose = 0;
Yold=zeros(size(Y));
for n=1:5
piter=1;
%result=ones(length(piter),4);

%for pi=1:length(piter)
pi=1
p=piter(pi);
[ii,jj]=find(Omega<=p);
YOmega=Y(Omega<=p);
[YOmega,ii,Jcol]=processInput(ii,jj,YOmega);
eps=zeros(d2,1);
for j=1:length(Jcol)-1
    ind = Jcol(j)+1:Jcol(j+1);
    if length(ind)<2
        continue
    end
    eps(j) = min(diff(YOmega(ind)));        
end


%% X = RMC_exact_fixed_margin(ii,jj,Jcol,YOmega,eps,d1,d2)
% min ||X||_* st DX_j<= -eps_j 

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

%% Initialize Variables
sv=20; 
t=1;%stepsize

Yrt=YOmega;
rinit=10;X.U=randn(d1,rinit);X.V=randn(d2,rinit);
XOmega=Amap(X);
spZ=ATmap(t*(Yrt-XOmega)); 

 %% Tracking
XOld=XOmega;YOld=Yrt; % for exit conditions



%% Continuation
mu0=1; 
result=[];   
mu=mu0;

for j=0:50
    mu=1^(-j)*mu0;
    fprintf('mu:%f\n',mu);
    ff=0.5*norm(XOmega-Yrt)^2;
    disp(ff)
    %% Iterations
    for iter=1:par.maxiter
        %% UPDATE                  
        %Xupdate    
        sv=SVT_LR_SP(t*mu,min(sv,par.maxrank));%SVT of X+spZ   
        %X.U=(sqrt(d1*d2/(sum(X.U.^2)*sum(X.V.^2)')))*X.U;
        XOmega=Amap(X);
        %Yupdate      
        if (rmc)
            Yrt=c_colMR_fixed_margin(XOmega',eps',Jcol'); Yrt=Yrt';
        end
        spZ=ATmap(t*(Yrt-XOmega));     

        %% EXIT CONDITIONS
        ch=norm(XOld-XOmega)/sqrt(length(XOmega));         
        ff=[0.5*norm(XOmega-Yrt)^2,mu*sum(sum(X.U.^2))];
        if par.verbose
            fprintf('\t iter:%d,sv:%d,ch1:%f, ch2:%0.2g,fval:%0.5g:%0.5g:%0.5g/%0.2g\n',...
                iter,sv,ch,norm(Yrt-YOld),sum(ff),ff(1),ff(2),norm(spZ,'fro'))            
        end  
        k=evalRanking(theta,X.U*X.V',f);k(3)=sqrt(k(3));%cmt            
        result=[result;[k',norm(spZ,'fro')]];%cmt
        XOld=XOmega;YOld=Yrt;
        
        if ((ch<par.tol) || (iter==par.maxiter))          
            break
        end                
    end
    if par.verbose
        Yrmc=X.U*X.V'+full(spZ);
        k=evalRanking(theta,Yrmc,f);k(3)=sqrt(k(3));
        k1=evalRanking(theta,X.U*X.V',f);k(3)=sqrt(k(3));
        fprintf('mu:%f, nucNorm:(%f,%f), normspZ:%f \n\t niter:%d ktau:(%f,%f), srho:(%f,%f), rmse:(%f,%f)\n',...
            mu,sum(svd(Yrmc)),sum(svd(X.U*X.V')),norm(spZ,'fro'),iter,k(1),k1(1),k(2),k1(2),k(3),k1(3));    
    end
    if (norm(spZ,'fro')<par.tol)
         break
    end    
end

%% 
Yrmc=X.U*X.V'+full(spZ);
norm(Yold-Yrmc)
Yold=Yrmc;
k=evalRanking(theta,Yrmc,f);k(3)=sqrt(k(3));
k1=evalRanking(theta,X.U*X.V',f);k(3)=sqrt(k(3));
fprintf('Size: %dX%d, rk:%d, p:%f, mu:(%f,%f), normspZ:%f \n\t iter:%d ktau:(%f,%f), srho:(%f,%f), rmse:(%f,%f)\n',...
        d1,d2,r,p,sum(svd(Yrmc)),sum(svd(X.U*X.V')),norm(spZ,'fro'),iter,k(1),k1(1),k(2),k1(2),k(3),k1(3));    
s=svd(Yrmc);    
end
%result(pi,1:3)=k; result(pi,4)=find(cumsum(s.^2)/sum(s.^2)>(1-1e-5),1); 
