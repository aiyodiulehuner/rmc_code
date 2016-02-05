%% X = RMC_exact_fixed_margin(ii,jj,Jcol,YOmega,eps,d1,d2)
% min ||X||_* st DX_j<= -eps_j 
function [Xest,spZest,stat]=RMC_exact_fixed_margin(ii,jj,Jcol,YOmega,eps,d1,d2,par,debug)
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
Yrt=YOmega;
rinit=10;
X.U=zeros(d1,rinit);X.V=zeros(d2,rinit); 
XOmega=Amap(X);
spZ=ATmap(Yrt-XOmega); 
Xold=XOmega;


%% Continuation
stat.result=[]; stat.iter=[]; stat.mu=[];
mu=par.mu0/par.continuation; 
res=0;
ch=0;
for j=0:50
    if res<par.tol && ch<par.tol
        mu=par.continuation*mu;    
        if ismember('mutarget', fieldnames(par))
            if mu<par.mutarget; mu=par.mutarget; end
        end
    else
        mu=sum(svd(X.U*X.V'+full(spZ))); par.continuation=1; 
    end
    if par.verbose; fprintf('mu:%f,%f\n',mu,0.5*norm(XOmega-Yrt)^2); end
    %% Iterations
    for iter=1:par.maxiter
        %% UPDATE                  
        sv=NNP_LR_SP(mu,min(sv,par.maxrank));%SVT/NNP of X+spZ   
        XOmega=Amap(X);
        Yrt=c_colMR_fixed_margin(XOmega',eps',Jcol); Yrt=Yrt';
        spZ=ATmap(Yrt-XOmega);     

        %% EXIT CONDITIONS
        res=norm(Yrt-XOmega)/sqrt(length(XOmega));   
        ch=norm(Xold-XOmega)/sqrt(length(XOmega));
        Xold=XOmega;
        if par.verbose
            fprintf('\titer:%d,sv:%d,res:%f/%0.2g,ch:%f,muY:%f\n',...
                iter,sv,res,ch,norm(spZ,'fro'),sum(svd(X.U*X.V'+full(spZ))))            
        end  
        if ~isempty(debug)
            theta=debug{1};f=debug{2};
            k=evalRanking(theta,X.U*X.V',f);k(3)=sqrt(k(3));
            stat.result=[stat.result;[k',norm(spZ,'fro')]];
        end

        if (res<par.tol && ch<par.tol)
            break
        end                
    end
    stat.iter=[stat.iter,iter];
    stat.mu=[stat.mu,mu];
    if ~isempty(debug)
        Yrmc=X.U*X.V'+full(spZ);
        k=evalRanking(theta,Yrmc,f);k(3)=sqrt(k(3));
        k1=evalRanking(theta,X.U*X.V',f);k(3)=sqrt(k(3));
        fprintf('mu:%f, nucNorm:(%f,%f), normspZ:%f \n\t niter:%d ktau:(%f,%f), srho:(%f,%f), rmse:(%f,%f)\n',...
            mu,sum(svd(Yrmc)),sum(svd(X.U*X.V')),norm(spZ,'fro'),iter,k(1),k1(1),k(2),k1(2),k(3),k1(3));    
    end
    if (par.continuation>=1)
         break
    end    
end
Xest=X;
spZest=spZ;
clear global
