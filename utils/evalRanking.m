function k=evalRanking(Xtrue,Yest,f)
% function to compute set of metrics over each column of Xtrue and Yest.
% f: cell array of function handles to the metrics

axis=2;
assert(all(size(Xtrue)==size(Yest)))
k=zeros(length(f),1);
for j=1:length(f)
    if ischar(f{j})
        f{j}=str2func(f{j});
    end
    for i=1:size(Xtrue,axis)
        k(j) = k(j) + f{j}(Xtrue(:,i),Yest(:,i));
    end
end
k=k/size(Xtrue,2);
end

function rho=spearman_rho(x,y)
    xrank=tiedrank(x,0);
    yrank=tiedrank(y,0);
    n=length(x);nconst=n*(n^2-1)/6;
    rho=1-(sum((xrank-yrank).^2)/nconst);
end

function tau=kendall_tau(x,y)
    xrank=tiedrank(x,0);
    yrank=tiedrank(y,0);
    n=length(x);nconst=n*(n-1)/2;
         
    tau = 0;
    for k = 1:n-1
        tau = tau + sum(sign(xrank(k)-xrank(k+1:n)).*sign(yrank(k)-yrank(k+1:n)));
    end
    tau = tau ./nconst;
    
end