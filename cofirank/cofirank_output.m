%sed 's/^/1 /' M.lsvm > M1.lsvm
%sed 's/^/1 /' U.lsvm > U1.lsvm
%ipython
%import numpy as np
%from sklearn.datasets import load_svmlight_file
%V,y=load_svmlight_file('U1.lsvm',100)
%U,y=load_svmlight_file('M1.lsvm',100)
%U=np.array(U.todense())
%V=np.array(V.todense())
%np.savetxt('U.csv', U, fmt='%f', delimiter=',', newline='\n')
%np.savetxt('V.csv', V, fmt='%f', delimiter=',', newline='\n')

%clear
%clc
addpath('utils/');
addpath('NNLS-0/solver');
addpath('NNLS-0/PROPACKmod/');

neurosynth=0;
movielens=1;
if neurosynth
K=10;
th=0.5;
probiter=0.2:0.2:0.8;
niter=1;
lambdaiter=5;
f={'spearman_rho', 'kendall_tau', 'NDCG', 'MSE', 'Precision'};
r=100;

resultCOFI=zeros(niter,length(probiter), length(lambdaiter), 3, length(f));
for i=1:niter
    for pi=1:probiter
        p=probiter(pi);        
        load(sprintf('../neurosynth_counts/folds/neurosynth_%d.mat',round(p*100)));
        fprintf('Size: %dX%d, p:%f, train:val:test::%d:%d:%d\n',d1,d2,p,...        
                 length(yy),length(yy_val),length(yy_test));
        for l=1:length(lambdaiter)            
            %parse files
            Ycofi.U=impoerdata(sprintf('cofirank_ns_%d/U.csv',round(p*100)));
            Ycofi.V=impoerdata(sprintf('cofirank_ns_%d/V.csv',round(p*100)));
            
            
            yest=Amap_MatComp(Ycofi,ii,Jcol);            
            k1=evalRanking(yy,yest,Jcol,f,K,th);            
            resultCOFI(i,pi,l,1,:)=k1;
            %validation            
            yest_val=Amap_MatComp(Ycofi,ii_val,Jcol_val);
            k2=evalRanking(yy_val,yest_val,Jcol_val,f,K,th);
            resultCOFI(i,pi,l,2,:)=k2;

            %testing
            yest_test=Amap_MatComp(Ycofi,ii_test,Jcol_test);
            k3=evalRanking(yy_test,yest_test,Jcol_test,f,K,th);                
            resultCOFI(i,pi,l,3,:)= k3;
            
            save('resultCOFI.mat','resultCOFI')
            save(sprintf('ycofi_%d_%d.mat',round(p*100),l),'Yrmc','t')
        end
    end
end
end
if movielens
cviter=5;
K=5;
th=3;
lambdaiter=5;
f={'spearman_rho', 'kendall_tau', 'NDCG', 'MSE', 'Precision'};

resultCOFI=zeros(cviter, length(lambdaiter), 3, length(f));
for cv=1:cviter
    load(sprintf('../data/ml-100k/folds/ml_%d.mat',cv));
    fprintf('Size: %dX%d, p:%f, train:val:test::%d:%d:%d\n',d1,d2,cv,...        
                length(yy),length(yy_val),length(yy_test));
    for l=1:length(lambdaiter)
        Ycofi.U=importdata(sprintf('results/cofirank/cofirank_ml_%d/U.csv',cv));
        Ycofi.V=importdata(sprintf('results/cofirank/cofirank_ml_%d/V.csv',cv));
        disp(Ycofi)
        
        %training
        yest=Amap_MatComp(Ycofi,ii,Jcol);            
        k1=evalRanking(yy,yest,Jcol,f,K,th);            
        resultCOFI(cv,l,1,:)=k1;

        %validation            
        yest_val=Amap_MatComp(Ycofi,ii_val,Jcol_val);
        k2=evalRanking(yy_val,yest_val,Jcol_val,f,K,th);
        resultCOFI(cv,l,2,:)=k2;

        %testing
        yest_test=Amap_MatComp(Ycofi,ii_test,Jcol_test);
        k3=evalRanking(yy_test,yest_test,Jcol_test,f,K,th);                
        resultCOFI(cv,l,3,:)= k3;
            
        save('results/ml_resultCOFI.mat','resultCOFI')
        save(sprintf('results/cofirank/ml_yrcofi_%d.mat',cv),'Ycofi')
    end
end
end
