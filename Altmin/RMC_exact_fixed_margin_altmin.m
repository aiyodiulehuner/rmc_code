%% X = RMC_exact_fixed_margin(ii,jj,Jcol,YOmega,eps,d1,d2)
% min ||X||_* st DX_j<= -eps_j
function [X,stat]=RMC_exact_fixed_margin_altmin(ii,jj,Jcol,YOmega,eps,d1,d2,par,debug)
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
Yrt=YOmega;
X.U=randn(d1,par.r);X.V=randn(d2,par.r);
XOmega=Amap(X);
Xold=XOmega;


stat.result=[]; stat.iter=[]; stat.mu=[];

if par.verbose; fprintf('mu:%f,%f\n',mu,0.5*norm(XOmega-Yrt)^2); end
%% Iterations
for iter=1:par.maxiter
    %% UPDATE
    X=AltMin_LR_SP(X,ATmap(Yrt));%SVT/NNP of X+spZ
    XOmega=Amap(X);
    Yrt=c_colMR_fixed_margin(XOmega',eps',Jcol); Yrt=Yrt';
    %% EXIT CONDITIONS
    res=norm(Yrt-XOmega)/sqrt(length(XOmega));
    ch=norm(Xold-XOmega)/sqrt(length(XOmega));
    Xold=XOmega;
    if par.verbose
        fprintf('\titer:%d,res:%f/%0.2g,ch:%f,muY:%f\n',...
            iter,res,ch,norm(spZ,'fro'),sum(svd(X.U*X.V')))
    end
    if ~isempty(debug)
        theta=debug{1};f=debug{2};
        k=evalRanking(theta,X.U*X.V',f);k(3)=sqrt(k(3));
        stat.result=[stat.result;[k',sum(svd(X.U*X.V'))]];
    end
    
    if (res<par.tol && ch<par.tol)
        break
    end
end
stat.iter=iter;

if ~isempty(debug)
    Yrmc=X.U*X.V';
    k=evalRanking(theta,Yrmc,f);k(3)=sqrt(k(3));
    k1=evalRanking(theta,X.U*X.V',f);k(3)=sqrt(k(3));
    fprintf('mu:%f, nucNorm:(%f,%f), normspZ:%f \n\t niter:%d ktau:(%f,%f), srho:(%f,%f), rmse:(%f,%f)\n',...
        mu,sum(svd(Yrmc)),sum(svd(X.U*X.V')),norm(spZ,'fro'),iter,k(1),k1(1),k(2),k1(2),k(3),k1(3));
end
