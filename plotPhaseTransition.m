%% Repeat
niter=5;
dim_iter={{100,100,5},{250,250,5}};
c=1:5;
giter = cell(0);
giter{1}=@(a) a;
for ci=1:length(c)   
    giter{ci+1}=@(a) 1.0./(1+exp(-c(ci)*a));
end
f={'spearman_rho', 'kendall_tau', 'NDCG'};
%%
muiter=[0.01,0.05,0.1,0.5,1,5,10,20];
probiter=0.05:0.05:0.95;  
%%
resultbest=zeros(niter, length(dim_iter), length(giter), length(probiter),length(f));
for ix=1:niter*length(dim_iter)
    [n,i]=ind2sub([niter,length(dim_iter)],ix);
    for gi=1:length(giter)
        for pi=1:length(probiter)
            temp=squeeze(result{ix}(gi,pi,:,:));
            [k,m]=max(temp);
            if length(unique(m))>1
                fprintf('Warning: argmax(%d,%d,%d)=',ix,gi,pi)
                disp(m)
            end
            resultbest(n,gi,i,pi,:)= k;
        end
    end
end

for gi=1:length(giter)
    h=figure('Name',func2str(giter{gi}));
    txt=cell(length(dim_iter),1);
    for i=1:length(dim_iter)
        [d1,d2,r]=dim_iter{i}{:};        
        for k=1:length(f)
            subplot(length(f),2,2*k-1);
            plot(probiter*(d1*d2)/((d1+d2)*(r-1)),...
                    mean(squeeze(resultbest(:,gi,i,:,k))),'-o'); hold on;
            xlabel('$\frac{|\Omega|}{2d(r-1)}$','Interpreter','latex')
            ylabel('evaluation')           
            subplot(length(f),2,2*k);
            plot(probiter*(d1*d2)/((d1+d2)*r*log(d1+d2)),...
                    mean(squeeze(resultbest(:,gi,i,:,k))),'-o'); hold on;
            xlabel('$\frac{|\Omega|}{2drlog(d)}$','Interpreter','latex')
            ylabel('evaluation')           
            txt{i}=num2str(dim_iter{i}{1});
        end
    end
    legend(txt)
end

for gi=1:length(giter)
    figure()
    for i=1:length(dim_iter)
        [d1,d2,r]=dim_iter{i}{:};
        for k=1:length(f)            
            plot(probiter*(d1*d2)/((d1+d2)*(r-1)),...
                squeeze(mean(resultbest(:,gi,i,:,k),1)),'-o'); hold on;
            xlabel('$\frac{|\Omega|}{2d(r-1)}$','Interpreter','latex')
            ylabel('evaluation')           
        end
    end
    legend({'Spearman Rho','Kendall Tau','NDCG'})
end
