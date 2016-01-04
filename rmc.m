function [Yest,iter,res]=rmc(ii,Jcol,jj,YOmega,d1,d2,mu,par) 

    %TODO: Line search, PAVPROBLEM CHECK
    %TODO: Continuation
    
    Amap  = @(X) Amap_MatComp(X,ii,Jcol);
    if (exist('mexspconvert')==3); 
        ATmap = @(y) mexspconvert(d1,d2,y,ii,Jcol); 
    else
        ATmap = @(y) spconvert([ii,jj,y; d1,d2,0]); 
    end
    
    % min_{X,Z,\epsilon} mu*\|X\|_*+0.5*\|P_Omega(X-Z)\|_2^2 s.t (Z,\epsilon)\in
    % \mathcal{X}_\Omega
    global X spZ
    
    sv=50;options.p0=rand(d1,1)-0.5;
    t=1;%stepsize
    
    % initialization
    rinit=10;
    X.U=zeros(d1,rinit);X.V=zeros(d2,rinit);
    spZ = ATmap(YOmega);
    [X.U,S,X.V]=lansvd('Axz','Atxz',d1,d2,rinit,'L',options);
    X.U=X.U*sqrt(S);
    X.V=X.V*sqrt(S);
    XOmega=Amap(X);
    
    %
    XOld=XOmega;
    Yrt=YOmega;
    eps=ones(d2,1);
    spZ=ATmap(t*(Yrt-XOmega)); 
    
    for iter=1:par.maxiter
        %Xupdate
        sv=SVT_LR_SP(t*mu,min(sv,par.maxrank),options);%SVT of X+spZ
        %Yupdate
        [Yrt,eps,stat]=RetargetColScores_margin(XOmega,eps,par.PAV_QP);
        if par.verbose
            fprintf('Retargeting exited with fval:%f, exitflag:%f\n',stat.fval,stat.exitflag)
        end
        options.p0=X.U(:,1);
        XOmega=Amap(X);        
        spZ=ATmap(t*(Yrt-XOmega));
        
        res=norm(XOld-XOmega);
        if res<par.tol
            if par.verbose
                fprintf('iter:%d,res:%2g,sv:%d,eps:%2g\n',iter,res,sv,min(eps))
            end
            break
        end
        
        XOld=XOmega;
        if par.verbose
            fprintf('iter:%d,res:%2g,sv:%d,eps:%2g\n',iter,res,sv,eps)
        end
    end
        
    Yest=X;

    clear global
    