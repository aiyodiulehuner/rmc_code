% min_{X,Z,\epsilon} mu*\|X\|_*+0.5*\|P_Omega(X-Z)\|_2^2 s.t (Z,\epsilon)\in
    % \mathcal{X}_\Omega
    

%TODO: Line search, PAVPROBLEM CHECK
%TODO: Continuation
    
function [Yest,iter,res,result]=rmc(ii,Jcol,jj,YOmega,d1,d2,mu,par,theta,f) 

    
    Amap  = @(X) Amap_MatComp(X,ii,Jcol);
    
    if (length(YOmega)/(d1*d2)>0.6)
        ATmap = @(y) full(sparse(d1,d2,y,ii,Jcol)); 
        
    else
        if (exist('mexspconvert')==3); 
            ATmap = @(y) mexspconvert(d1,d2,y,ii,Jcol); 
        else
            ATmap = @(y) sparse([ii,jj,y; d1,d2,0]); 
        end
    end
    
    
    global X spZ
    
    sv=50;options.p0=rand(d1,1)-0.5;
    t=1;%stepsize
    
    % initialization
    rinit=10;
   
    eps=zeros(d2,1);
    for j=1:length(Jcol)-1
        ind = Jcol(j)+1:Jcol(j+1);
        eps(j) = min(diff(YOmega(ind)));
    end
    assert(all(eps>=0));
    scale=sum(eps);
    Yrt = YOmega/scale;    
    eps=eps/scale; 
    mu=mu/scale;
    
    X.U=zeros(d1,rinit);X.V=zeros(d2,rinit);
    spZ = ATmap(Yrt);  
    [X.U,S,X.V]=lansvd('Axz','Atxz',d1,d2,rinit,'L',options);
    X.U=X.U*sqrt(S);
    X.V=X.V*sqrt(S);
    XOmega=Amap(X);
    
    XOld=XOmega;
    YOld=Yrt;
    XXOld=X;
    epsOld=eps;
    
    spZ=ATmap(t*(Yrt-XOmega)); 
    fold=0.5*norm(XOmega-Yrt)^2+mu*sum(sum(X.U.^2));
    disp(fold)
    result=[];
    for iter=1:par.maxiter
        %Xupdate
        sv=SVT_LR_SP(t*mu,min(sv,par.maxrank));%SVT of X+spZ
        %Yupdate
        [Yrt,eps,stat]=RetargetColScores_margin(XOmega,eps,par.PAV_QP);
        if par.verbose && abs(stat.exitflag-1)>1e-10
            fprintf('Retargeting exited with fval:%f, exitflag:%f\n',stat.fval,stat.exitflag)
        end        
        XOmega=Amap(X);        
        spZ=ATmap(t*(Yrt-XOmega));
        
        ch=norm(XOld-XOmega)/sqrt(length(XOmega));     
        res=norm(XOmega-Yrt)/sqrt(length(XOmega));
        if ((ch<par.tol) || (res<par.tol))
            if par.verbose
                fprintf('iter:%d,res:%2g,sv:%d,eps:%2g, ch1:%f, ch2:%0.2g,fval:%0.2g\n',...
                    iter,res,sv,min(eps),ch,norm(Yrt-YOld)/sqrt(length(XOmega)),0.5*norm(XOmega-Yrt)^2+mu*sum(sum(X.U.^2)))

            end
            break
        end
        
        ff=0.5*norm(XOmega-Yrt)^2+mu*sum(sum(X.U.^2));
        if par.verbose
            fprintf('iter:%d,res:%2g,sv:%d,eps:%2g, ch1:%f, ch2:%0.2g,fval:%0.5g\n',...
                    iter,res,sv,min(eps),ch,norm(Yrt-YOld),ff)
            k=evalRanking(theta,X.U*X.V',f);k(3)=sqrt(k(3));
            result=[result;[k',nnz(eps)]];
        end
        if (fold<ff)
           fprintf('Warning:non-descent step')
            break
        end    
        
        XOld=XOmega;
        YOld=Yrt;
        fold=ff;
        XXOld=X;
        epsOld=eps;
    
    end
        
    Yest=X;

    clear global
 