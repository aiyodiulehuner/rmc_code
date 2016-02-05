function resultSMC=runSMCSimulation(d1,d2,r,giter,f)
    
   %parameters
    par.tol     = 1e-5;
    par.maxiter = 1000;
    par.maxrank = 100;%min([d1,d2,500]);
    par.verbose = 0;
    par.nnp=1;

    muiter=[0.05,0.1,0.5,1,5,10];
    probiter=0.05:0.05:0.95;    
    %probiter=[0.4];
    resultSMC=zeros(length(giter), length(probiter), length(muiter), length(f));
    U=randn(d1,r);
    V=randn(d2,r);
    theta=U*V'; 
    theta=sqrt(d1*d2)*theta/norm(theta,'fro');
    Omega=rand(size(theta));
    
    mu0=sum(svd(theta));
    
    for gi=1:length(giter)
        g=giter{gi};
        Y=theta;
        for j=1:d2
            Y(:,j)=g(theta(:,j));
        end

        for pi=1:length(probiter)
            p=probiter(pi);
            
            [ii,jj]=find(Omega<=p);
            YOmega=Y(Omega<=p);
            [YOmega,ii,Jcol]=processInput(ii,jj,YOmega);
            fprintf('Size: %dX%d, rk:%d, p:%f, gi:%d\n',d1,d2,r,p,gi);
                
            for m=1:length(muiter)
                mu=mu0*muiter(m);
                [Ysmc,iter,res]=smc(ii,Jcol,jj,YOmega,d1,d2,mu,par);
                k=evalRanking(theta,Ysmc.U*Ysmc.V',f);k(3)=sqrt(k(3));
                fprintf('\t mu:%f. iter:%d, res:%f, ||X||_*:%f, ktau:%f, srho:%f, ndcg:%f\n',...
                    mu, iter,res,sum(sum(Ysmc.U.^2)),k(1),k(2),k(3));
                resultSMC(gi,pi,m,:)= k;                
            end

        end          
    end    
