function k=evalRanking(ytrue,yest,Jcol,f,K,th)
% function to compute set of metrics over each column of Xtrue and Yest.
% f: cell array of function handles to the metrics
d2=length(Jcol)-1;
assert(all(size(ytrue)==size(yest)))
k=zeros(length(f),1);
nd=0;
for i=1:d2
    ind=Jcol(i)+1:Jcol(i+1);
    if (~isempty(ind)) && (length(unique(yest(ind)))>1) && (length(unique(ytrue(ind)))>1)
        nd=nd+1;
        for j=1:length(f)
            if ischar(f{j})
                f{j}=str2func(f{j});
            end    
            k(j) = k(j) + f{j}(ytrue(ind),yest(ind),K,th);
        end
    end
end
k=k/nd;
end

function prec_k=Precision(x,y,K,th)
   xth=(x>th);
   [yy,ii]=sort(y,'descend');
   l=xth(ii);
   r=(1:length(yy))';
   all_prec=cumsum(l)./r;
   if K<length(x)
      prec_k=all_prec(K);
   else
      prec_k=all_prec(end);
   end
end


function mse=MSE(x,y,K,th)
   mse=norm(x-y)^2/length(x);
end

%y is score and x is ground truth
function rho=spearman_rho(x,y,K,th)
    %xrank=tiedrank(x,0);
    %yrank=tiedrank(y,0);
    %n=length(x);nconst=n*(n^2-1)/6;
    %rho=1-(sum((xrank-yrank).^2)/nconst);
    if size(x,2)>1
        x=x';
    end
    if size(y,2)>1
        y=y';
    end
    if (sum(abs(x))==0)
        disp('zero x')
    end    
    rho = corr(x,y,'type','Spearman');
end

function tau=kendall_tau(x,y,K,th)
    %xrank=tiedrank(x,0);
    %yrank=tiedrank(y,0);
    %n=length(x);nconst=n*(n-1)/2;
         
    %tau = 0;
    %for k = 1:n-1
    %    tau = tau + sum(sign(xrank(k)-xrank(k+1:n)).*sign(yrank(k)-yrank(k+1:n)));
    %end
    %tau = tau ./nconst;   
    if size(x,2)>1
        x=x';
    end
    if size(y,2)>1
        y=y';
    end
    tau = corr(x,y,'type','Kendall');
end

function ndcg_k=NDCG(x,y,K,th)
%x.y need to be columns    
    if size(x,2)>1
        x=x';
    end
    if size(y,2)>1
        y=y';
    end

    [yy,ii]=sort(y,'descend');
    r=(1:length(yy))';%r=1 for highest y, 2 for second highest, ...ri=rank(i) in y
    l=x(ii);%l=l/max(l); % true relevance of item in rth position 
    log2r = log2(r+2.0);
    exp2l = (2.0.^l) - 1.0;
    [~,indg]  = sort(l,'descend');
        
    dcg_k  = exp2l./log2r; % numerator dcg
    dcg_k = cumsum(dcg_k);
    idcg_k = exp2l(indg)./log2r;  % ideal dcg
    idcg_k = cumsum(idcg_k);
    idcg_k(idcg_k==0)=1.0;
    ndcg = dcg_k(end)/idcg_k(end);
    all_ndcg = dcg_k./idcg_k;
        
    if K > length(x)
        ndcg_k =ndcg;            
    else
        ndcg_k = all_ndcg(K);
    end
end

