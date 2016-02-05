%% Repeat
niter=1;
dim_iter={{100,100,5}};
c=2:2:8;
giter = cell(0);
giter{1}=@(a) a;
for ci=1:length(c)   
    giter{ci+1}=@(a) 1.0./(1+exp(-c(ci)*a));
end
f={'Spearman \rho', 'Kendall \tau', 'NDCG'};
%%
muiter=[1e-4,1e-3,0.01,0.05];
probiter=0.1:0.1:0.9;  
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
            %k(3)=k(3)^2;
            resultbest(n,i,gi,pi,:)= k;
        end
    end
end
figure()

for gi=1:length(giter)    
    for i=1:length(dim_iter)
        [d1,d2,r]=dim_iter{i}{:};        
        for k=1:length(f)
            subplot(1,length(f),k);
            plot(probiter*(d1*d2)/((d1+d2)*r*log(d1+d2)),...
                    (squeeze(resultbest(:,i,gi,:,k))),'-o','LineWidth',2); hold on;
            xlabel('$\frac{|\Omega|}{2drlog(d)}$','Interpreter','latex','FontSize',16)
            ylabel(f(k),'FontSize',16)           
            txt{gi}=[num2str(dim_iter{i}{1}),'-',func2str(giter{gi})];
        end
    end    
end
legend({'\Theta','1/(1+exp(-2\Theta))','1/(1+exp(-4\Theta))','1/(1+exp(-6\Theta))','1/(1+exp(-8\Theta))'})

% for gi=1:length(giter)
%     figure()
%     for i=1:length(dim_iter)
%         [d1,d2,r]=dim_iter{i}{:};
%         for k=1:length(f)            
%             plot(probiter*(d1*d2)/((d1+d2)*(r-1)),...
%                 squeeze(mean(resultbest(:,gi,i,:,k),1)),'-o'); hold on;
%             xlabel('$\frac{|\Omega|}{2d(r-1)}$','Interpreter','latex')
%             ylabel('evaluation')           
%         end
%     end
%     legend({'Spearman Rho','Kendall Tau','NDCG'})
% end
