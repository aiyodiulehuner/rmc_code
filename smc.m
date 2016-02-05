%% STANDARD MATRIX COMPLETION
% \min_X \|P_\Omega(X-Y)\|^2 s.t. \|X\|_*\le mu <- NNP_LR_SP
% \min_X 0.5*\|P_\Omega(X-Y)\|^2 + \mu*\|X\|_*

function [Yest,iter,res]=smc(ii,Jcol,jj,YOmega,d1,d2,mu,par) 
    Amap  = @(X) Amap_MatComp(X,ii,Jcol);
    if (exist('mexspconvert')==3); 
        ATmap = @(y) mexspconvert(d1,d2,y,ii,Jcol); 
    else
        ATmap = @(y) spconvert([ii,jj,y; d1,d2,0]); 
    end    
    
    global X spZ
    sv=50;
    
    X.U=zeros(d1,10);X.V=zeros(d2,10);
    XOmega=zeros(size(YOmega));
    spZ=ATmap(YOmega-XOmega);                    
    XOld=XOmega;
    for iter=1:par.maxiter
        %disp(sv)
        %disp(min(sv,par.maxrank))
        if (ismember('nnp',fieldnames(par)) && par.nnp==1)
            sv=NNP_LR_SP(mu,min(sv,par.maxrank));
        else
            sv=SVT_LR_SP(mu,min(sv,par.maxrank));
        end
        
        XOmega=Amap(X);
        res=norm(XOmega-YOmega);
        ch=norm(XOld-XOmega);
        if ch<par.tol || res<par.tol 
            break
        end
        
        XOld=XOmega;
        spZ=ATmap(YOmega-XOmega);
        
        if par.verbose
            fprintf('iter:%d,res:%2g,ch=%f,sv:%d, mu:%f/%f\n',iter,res,ch,sv,mu,sum(sum(X.U.^2)))
        end
    end
        
    Yest=X;
    clear global
    