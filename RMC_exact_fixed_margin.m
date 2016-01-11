%% X = RMC_exact_fixed_margin(ii,jj,Jcol,YOmega,eps,d1,d2)
% min ||X||_* st DX_j<= -eps_j 
function RMC_exact_fixed_margin(ii,jj,Jcol,YOmega,eps,d1,d2,debug)
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
X.U=randn(d1,rinit);X.V=randn(d2,rinit); 
XOmega=Amap(X);
spZ=ATmap(Yrt-XOmega); 


%% Continuation
continuation=0.5;
result=[];   
mu=mu0/continuation; 
ch=0;
for j=0:50
    if ch<par.tol; mu=continuation*mu;    
    else mu=1.1*mu; continuation=1.1*continuation; end
    fprintf('mu:%f,%f\n',mu,0.5*norm(XOmega-Yrt)^2);
    %% Iterations
    for iter=1:par.maxiter
        %% UPDATE                  
        sv=NNP_LR_SP(mu,min(sv,par.maxrank));%SVT/NNP of X+spZ   
        XOmega=Amap(X);
        Yrt=c_colMR_fixed_margin(XOmega',eps',Jcol'); Yrt=Yrt';
        spZ=ATmap(Yrt-XOmega);     

        %% EXIT CONDITIONS
        ch=norm(Yrt-XOmega)/sqrt(length(XOmega));         
        if par.verbose
            fprintf('\titer:%d,sv:%d,ch:%f/%0.2g\n',iter,sv,ch,norm(spZ,'fro'))            
        end  
        k=evalRanking(theta,X.U*X.V',f);k(3)=sqrt(k(3));
        result=[result;[k',norm(spZ,'fro')]];

        if (ch<par.tol)
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
    if (continuation>=1)
         break
    end    
end