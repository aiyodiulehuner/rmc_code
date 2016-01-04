%% STANDARD MATRIX COMPLETION
% \min_X \|P_\Omega(X-Y)\|^2 s.t. \|X\|_*\le mu <- NNP_LR_SP
% \min_X 0.5*\mu*\|P_\Omega(X-Y)\|^2 + \|X\|_*

function [Yest,iter,res]=smc(ii,Jcol,jj,YOmega,d1,d2,mu,par) 
    
    Amap  = @(X) Amap_MatComp(X,ii,Jcol);
    if (exist('mexspconvert')==3); 
        ATmap = @(y) mexspconvert(d1,d2,y,ii,Jcol); 
    else
        ATmap = @(y) spconvert([ii,jj,y; d1,d2,0]); 
    end
    
    
    global X spZ
    sv=50;options.p0=rand(d1,1)-0.5;
    
    X.U=zeros(d1,10);X.V=zeros(d2,10);
    XOmega=zeros(size(YOmega));
    XOld=XOmega;
    spZ=ATmap(YOmega-XOmega);                    
    
    for iter=1:par.maxiter
        sv=NNP_LR_SP(mu,min(sv,par.maxrank),options);
        options.p0=X.U(:,1);
        XOmega=Amap(X);
        res=norm(XOld-XOmega);
        if res<par.tol
            break
        end
        XOld=XOmega;
        spZ=ATmap(YOmega-XOmega);
        if par.verbose
            fprintf('iter:%d,res:%2g,sv:%d\n',iter,res,sv)
        end
    end
        
    Yest=X;

    clear global
    