%% X = RMC_exact_fixed_margin(ii,jj,Jcol,YOmega,eps,d1,d2)
% min ||X||_* st DX_j<= -eps_j 
%[Xest,spZest,stat]
function [Yest,iter,res]=rmc_fixed_margin(ii,Jcol,jj,YOmega,d1,d2,mu,par)

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

%% Initialize Variables
sv=20; 

global X spZ
rinit=10;
eps=zeros(d2,1);
for j=1:length(Jcol)-1
    ind = Jcol(j)+1:Jcol(j+1);
    eps(j) = max(1e-10,min(diff(YOmega(ind))));
end
%eps=max(eps,1e-10)*ones(d2,1);
Yrt=YOmega;

X.U=zeros(d1,rinit);X.V=zeros(d2,rinit); 
XOmega=Amap(X);
spZ=ATmap((Yrt-XOmega)/2);
Xold=XOmega;
par.continuation=0.5;mu0=222/par.continuation;
ch=0; res=0; mu=mu0;
for j=0:0
    if res<par.tol && ch<par.tol
        mu=par.continuation*mu;    
        if ismember('mutarget', fieldnames(par))
            if mu<par.mutarget; mu=par.mutarget; end
        end
    else
        mu=sum(svd(X.U*X.V'+full(spZ))); par.continuation=1; 
    end
    
for iter=1:par.maxiter
    %% UPDATE                  
    sv=NNP_LR_SP(mu,min(sv,par.maxrank));
    XOmega=Amap(X);
    Yrt=c_colMR_fixed_margin(((Yrt+XOmega)/2)',eps',Jcol); Yrt=Yrt';  

    spZ=ATmap((Yrt-XOmega)/2);     

    %% EXIT CONDITIONS
    res=norm(Yrt-XOmega)/sqrt(length(XOmega));   
    ch=norm(Xold-XOmega)/sqrt(length(XOmega));
    Xold=XOmega;
    if par.verbose
        fprintf('\titer:%d,sv:%d,res:%f/%0.2g,ch:%f,muY:%f\n',...
            iter,sv,res,norm(spZ,'fro'),ch,sum(svd(X.U*X.V'+full(spZ))))            
    end  

    if (res<par.tol || ch<par.tol)
        break
    end                
end
k=evalRanking(theta,X.U*X.V',f)
if (par.continuation>=1)
     break
end  
end
Yest=X;

clear global
