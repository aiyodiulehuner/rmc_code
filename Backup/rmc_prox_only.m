% min_{X,Z,\epsilon} mu*\|X\|_*+0.5*\|P_Omega(X-Z)\|_2^2 s.t (Z,\epsilon)\in
    % \mathcal{X}_\Omega
    

%TODO: Line search, PAVPROBLEM CHECK
%TODO: Continuation
    
function [Yest,iter,res,result]=rmc(ii,Jcol,jj,YOmega,d1,d2,mu,par,theta,f) 
    acceleration=0;
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
    
    sv=20;options.p0=rand(d1,1)-0.5;
    t=1;%stepsize
    
    %% initialization X0,Y0,eps0
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
    
    %% Assign Variables
    X.U=zeros(d1,rinit);X.V=zeros(d2,rinit);
    spZ = ATmap(Yrt);  
    [X.U,S,X.V]=lansvd('Axz','Atxz',d1,d2,rinit,'L',options);
    X.U=X.U*sqrt(S);
    X.V=X.V*sqrt(S);
    XOmega=Amap(X);
    spZ=ATmap(t*(Yrt-XOmega));     
    %----------[X,Y,eps]=[X,Y,eps]^(k),XOmega_=P_O(X^(k))   
        
    X_=X;
    Yrt_=Yrt;
    eps_=eps;
    %----------[X_,Y_,eps_]=[X,Y,eps]^(k-1)=[X,Y,eps]^(k)
    
    %% Tracking
    XOld=XOmega;YOld=Yrt; % for exit conditions
    ff=[0.5*norm(XOmega-Yrt)^2,mu*sum(sum(X.U.^2))];
    disp(ff)
    result=[];
    
    %% Iterations
    for iter=1:par.maxiter
        %% UPDATE
        %----------[X,Y,eps]=[Xa,Ya,eps_a]^(k+1),spZ=(Ya^(k+1)-P_O(Xa^(k+1))),XOmega=P_O(Xa^(k+1))              
        %Yupdate              
        [Yrt,eps,stat]=RetargetColScores_margin(XOmega,eps,par.PAV_QP);                
        if par.verbose && abs(stat.exitflag-1)>1e-10
            fprintf('Warning: Retargeting exited with fval:%f, exitflag:%f\n',stat.fval,stat.exitflag)
        end 
        %Xupdate
        spZ=ATmap(t*(Yrt-XOmega));     
        sv=SVT_LR_SP(t*mu,min(sv,par.maxrank));%SVT of X+spZ  
        XOmega=Amap(X);
        %----------[X,Y,eps]=[X,Y,eps]^(k+1), XOmega=P_O(X^(k+1))
        
        
        %% EXIT CONDITIONS
        ch=norm(XOld-XOmega)/sqrt(length(XOmega));     
        res=norm(XOmega-Yrt)/sqrt(length(XOmega));
        ff=[0.5*norm(XOmega-Yrt)^2,mu*sum(sum(X.U.^2))];
        if par.verbose
            fprintf('iter:%d,res:%2g,sv:%d,eps:%2g, ch1:%f, ch2:%0.2g,fval:%0.5g:%0.5g:%0.5g\n',...
                iter,res,sv,min(eps),ch,norm(Yrt-YOld),sum(ff),ff(1),ff(2))
            k=evalRanking(theta,X.U*X.V',f);k(3)=sqrt(k(3));%cmt
            result=[result;[k',min(eps)]];%cmt
        end     
        XOld=XOmega;YOld=Yrt;
        if ((ch<par.tol) || (res<par.tol) || (iter==par.maxiter))          
            break
        end                
        
        
        %% ACCELERATION (1+w)*[X,Y,eps]^(k+1)-w*[X,Y,eps]^(k)   
        if (acceleration)
            w=2/(iter+1); 
            w=0;       
            X_t.U=[(1+w)*X.U,-w*X_.U];X_t.V=[X.V,X_.V];
            eps_t=(1+w)*eps-w*eps_;
            Yrt_t=(1+w)*Yrt-w*Yrt_;
            %---------[X_t,Y_t,eps_t]=[Xa,Ya,eps_a]^(k+2)

            %% SETUP FOR NEXT ITERATION 
            X_=X;
            Yrt_=Yrt;
            eps_=eps;
            %--------[X_,Y_,eps_]=[X,Y,eps]^(k+1) --for (k+3) update

            X=X_t;
            Yrt=Yrt_t;
            eps=eps_t;

            XOmega=Amap(X);         
            spZ=ATmap(t*(Yrt-XOmega));     
            %---------[X,Y,eps]=[Xa,Ya,eps_a]^(k+2),spZ=(Ya^(k+2)-P_O(Xa^(k+2))),XOmega=P_O(Xa^(k+1))
        end
    end
        
    Yest=X;

    clear global
 