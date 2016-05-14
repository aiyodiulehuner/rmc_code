function k=evalRanking(Xtrue,Yest,f)
% function to compute set of metrics over each column of Xtrue and Yest.
% f: cell array of function handles to the metrics

axis=2;
assert(all(size(Xtrue)==size(Yest)))
k=zeros(length(f),1);
for i=1:size(Xtrue,axis)
    for j=1:length(f)
        if ischar(f{j})
            f{j}=str2func(f{j});
        end    
        k(j) = k(j) + f{j}(Xtrue(:,i),Yest(:,i));
    end
end
k=k/size(Xtrue,2);
end

%y is score and x is ground truth
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

function ndcg_k=NDCG(x,y)
%x.y need to be columns
    if size(x,2)>1
        x=x';
    end
    if size(y,2)>1
        y=y';
    end
    k=50;
    [xx,ii]=sort(x,'descend');
    r=(1:length(xx))';%r=1 for highest x, 2 for second highest, ...ri=rank(i) in x
    l=y(ii);%l=l/max(l);
    log2r = log2(r+1.0);
    exp2l = (2.0.^l) - 1.0;
    [~,indg]  = sort(l,'descend');
        
    dcg_k  = exp2l./log2r; % numerator dcg
    dcg_k = cumsum(dcg_k);
    idcg_k = exp2l(indg)./log2r;  % ideal dcg
    idcg_k = cumsum(idcg_k);
    idcg_k(idcg_k==0)=1.0;
    ndcg = dcg_k(end)/idcg_k(end);
    all_ndcg = dcg_k./idcg_k;
        
    if k > length(x)
        ndcg_k =ndcg;            
    else
        ndcg_k = all_ndcg(k);
    end
end