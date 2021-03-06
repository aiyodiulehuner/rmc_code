%% X = RMC_exact_fixed_margin(ii,jj,Jcol,YOmega,eps,d1,d2)
% min ||X||_* st DX_j<= -eps_j 
% [Xest,spZest,stat]
function [Yest,Yrt,iter,res,ii]=rmc_fixed_margin(ii,Jcol,jj,YOmega,d1,d2,mu0,par,Xinit,Yinit)
Amap  = @(X,ii) Amap_MatComp(X,ii,Jcol); 
if (length(YOmega)/(d1*d2)>0.6)
    ATmap = @(y,ii) full(sparse(ii,jj,y, d1,d2));
else
    if (exist('mexspconvert')==3); 
        ATmap = @(y,ii) mexspconvert(d1,d2,y,ii,Jcol); 
    else
        ATmap = @(y,ii) sparse(ii,jj,y, d1,d2);
    end
end

%% Initialize Variables
sv=par.maxrank;

global X spZ
rinit=10;

n=length(YOmega);
%compute epsilon
eps0=1e-3;
eps=zeros(n,1);
blk={};

for j=1:length(Jcol)-1
    ind = Jcol(j)+1:Jcol(j+1);
    Yj=diff(YOmega(ind));%diff(y)=y(i)-y(i+1)

    eps_temp=eps0*(Yj>=1e-5);
    eps(ind)=[0;cumsum(eps_temp)];
    eps(ind)=eps(ind)-max(eps(ind))/2.0;
    %create blks
    f= find(eps_temp(1:end)==0);
    if ~isempty(f)
        sid=f(1);
        for i=2:length(f)
            if f(i)==f(i-1)+1
                continue
            else
                eid=f(i-1)+1;
                if (eid-sid)>1
                    blk{length(blk)+1}=[ind(sid),ind(eid)];
                end
                sid=f(i);
            end
        end
        i=length(f);
        eid=f(i)+1;
        if (eid-sid)>1
           blk{length(blk)+1}=[ind(sid),ind(eid)];
        end
    end
end
fprintf('len(blk):%d,max(eps):%d\n',length(blk),max(eps))
eps=eps/max(eps);

Yrt=Yinit;%Omega;
X.U=Xinit.U;X.V=Xinit.V;
XOmega=Amap(X,ii);
eta=1;
spZ=ATmap(eta*(Yrt-XOmega),ii);
Xold=XOmega;
Yold=Yrt;
if mu0<0 
    continuation_steps=4;
else
    continuation_steps=1;
end
par.continuation=0.5;mu0=mu0/((par.continuation)^continuation_steps);
res=0; mu=mu0;
idx=1:length(Yrt);
for j=1:continuation_steps
    mu=par.continuation*mu;

    for iter=1:par.maxiter
        %% UPDATE 
        if par.nnp
            sv=NNP_LR_SP(mu,sv,par);

            XOmega=Amap(X,ii);
            Yrt_temp=(1-eta)*Yrt+eta*XOmega;  
            if ~isempty(blk)
                [Yrt_temp,ii,idx]=block_sort(Yrt_temp,ii,blk);
            end
            Yrt=c_colMR_fixed_margin(Yrt_temp',eps',Jcol'); Yrt=Yrt';
            spZ=ATmap(eta*(Yrt-XOmega),ii);
	    
	    chY=norm(Yrt-Yold(idx))^2/n;
	    fprintf(' Ych:%f\n',chY)
            chX=norm(XOmega-Xold)^2/n;
            fprintf('\t\t NNP: sv:%d, mu:%f, Xch:%f\n',sv,mu,chX)
            ch=max(chX,chY);
            
            Xold=XOmega;
            Yold=Yrt;
        else
            sv=SVT_LR_SP(mu,sv,par); 
            fprintf('\n\t\t SVT: sv:%d,muX:%f',sv,sum(svd(X.U*X.V'))); 
 
            ch=norm(Amap(X,ii)-Xold)^2/n;
            Xold=Amap(X,ii);
            
            Yrt_temp=(Yrt+XOmega)/2; 
            if ~isempty(blk)
                [Yrt_temp,ii]=block_sort(Yrt_temp,ii,blk);
            end
            Yrt=c_colMR_fixed_margin(Yrt_temp',eps',Jcol'); Yrt=Yrt';
            
            XOmega=Amap(X,ii); 
            spZ=ATmap((Yrt-XOmega)/2,ii);
        end 
        %% EXIT CONDITIONS
        res=norm(Yrt-XOmega);
        %ch=norm(Xold-XOmega)/sqrt(length(XOmega));
        %Xold=XOmega;
        if par.verbose
            fprintf('\titer:%d,sv:%d,res:%f/%f,ch:%f,muY:%f\n',...
                iter,sv,res,norm(spZ,'fro')/eta,ch,sum(sum(X.U.^2)))
        end
        if (res<par.tol || ch<par.tol^2)
            break
        end
    end
end
Yest=X;

clear global
