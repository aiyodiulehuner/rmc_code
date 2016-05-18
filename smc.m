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
    sv=1000;
    
    X.U=zeros(d1,10);X.V=zeros(d2,10);
    XOmega=zeros(size(YOmega));
    spZ=ATmap(YOmega-XOmega);                    
    XOld=XOmega;
    for iter=1:par.maxiter
        %disp(sv)
        %disp(min(sv,par.maxrank))
        if (ismember('nnp',fieldnames(par)) && par.nnp==1)           
            sv=NNP_LR_SP(mu,sv,par);
            fprintf('\t\tNNP: sv:%d, mu:%f\n',sv,mu)
        else            
            sv=SVT_LR_SP(mu,sv,par);
            fprintf('\t\t SVT: sv:%d,muX:%f\n',sv,sum(svd(X.U*X.V')));   
        end
        
        XOmega=Amap(X);
        res=norm(XOmega-YOmega);
        ch=norm(XOld-XOmega)/length(XOmega);
        if par.verbose
            fprintf('\titer:%d,sv:%d,res:%f/%f,ch:%f,muY:%f\n',...
                iter,sv,res,2.0*norm(spZ,'fro'),ch,sum(sum(X.U.^2)))              
        end
        if ch<par.tol || res<par.tol 
            break
        end
        
        XOld=XOmega;
        spZ=ATmap(YOmega-XOmega);                
    end
        
    Yest=X;
    clear global
    