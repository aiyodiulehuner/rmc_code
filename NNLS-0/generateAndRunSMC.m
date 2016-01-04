function generateAndRunSMC(d,r,prob)
    U = randn(d(i),r(i));
    V = randn(d(i),r(i));
    Omega = rand(d,d);
    theta = U*V';
    
    
    par.tol     = 1e-4;
    par.verbose = 0;
    par.truncation = 0; 
    par.continuation = 0;
    par.maxiter = 1000;
    par.maxrank = 100;
    

    
    for p = prob
        [ii,jj] = find(Omega<=p);
        YOmega=theta(Omega<=p);
        n=length(ii);
        Jcol=compJcol(jj);
        for j=1:d
            Yj=YOmega(Jcol(j)+1:Jcol(j+1));
            [~,IX]=sort(Yj,'ascend');
            a=zeros(length(IX),1); a(IX)=1:length(IX);
            YOmega(Jcol(j)+1:Jcol(j+1))=a;
        end
    end
end
        
        
        lamda=(d*d/n)*[1e-6,1e-5,1e-4,0.001,0.01,0.1,1,10,100,1000,1e4,1e5];
        param_iter=(1:12);
        
        Amap  = @(X) Amap_MatComp(X,ii,Jcol);
        if (exist('mexspconvert')==3); 
            ATmap = @(y) mexspconvert(d1,d2,y,ii,Jcol); 
        else
            ATmap = @(y) spconvert([ii,jj,y; d1,d2,0]); 
        end
        
        for pp=param_iter    
            param=lamda(pp);
            [X, iter, time, sd, hist] = APGL(d,d,'NNLS',Amap,ATmap,YOmega,1.0/param,0,par);
        
        
    if strcmp(g,'argsort')
        for j=1:d2
            Yj=YOmega(Jcol(j)+1:Jcol(j+1));
            [y,IX]=sort(Yj,'ascend');
            a=zeros(length(IX),1); a(IX)=1:length(IX);
            YOmega(Jcol(j)+1:Jcol(j+1))=a;