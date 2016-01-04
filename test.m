
addpath('utils/');
addpath('NNLS-0/solver');
addpath('NNLS-0/PROPACKmod/');
close all;
clc;clear;clear global;

rng(42)
funcSimulation=@runRMCSimulation;
outfile='resultRMC_test';
niter=1;
dim_iter={{100,100,5}};

giter = cell(0);
giter{1}=@(a) a;
f={'spearman_rho', 'kendall_tau', @(x,y)norm(x-y)^2/length(x-y)};

result=cell(niter*length(dim_iter),1);

parfor ix=1:niter*length(dim_iter)
    [n,i]=ind2sub([niter,length(dim_iter)],ix);
    [d1,d2,r]=dim_iter{i}{:};
    result{ix}=funcSimulation(d1,d2,r,giter,f);
end
save(outfile,'result')


fprintf('DONE COMPUTING\n')
load resultSMC_test
niter=30;

dim_iter={{100,100,5}};

giter = cell(0);
giter{1}=@(a) a;
f={'spearman_rho', 'kendall_tau', @(x,y)norm(x-y)^2/length(x-y)};


muiter=[0.01,0.05,0.1,0.2,0.5,0.7,1,5,10,20];
probiter=0.05:0.05:0.95;  
resultbest=zeros(niter, length(dim_iter), length(giter), length(probiter),length(f));
for ix=1:niter*length(dim_iter)
    [n,i]=ind2sub([niter,length(dim_iter)],ix);
    for gi=1:length(giter)
        for pi=1:length(probiter)
            temp=squeeze(result{ix}(gi,pi,:,:));
            temp(:,end)=-temp(:,end);
            [k,m]=max(temp);k(end)=-k(end);
            if length(unique(m))>1
                fprintf('Warning: argmax(%d,%d,%d)=',ix,gi,pi)
                disp(m)
            end
            resultbest(n,gi,i,pi,:)= k;
        end
    end
end
for gi=1:2:length(giter)
    figure()
    for i=1:length(dim_iter)
        [d1,d2,r]=dim_iter{i}{:};        
        for k=1:length(f)
            subplot(1,2,1);
            plot(probiter*(d1*d2)/((d1+d2)*(r-1)),...
                    mean(squeeze(resultbest(:,gi,i,:,k))),'-o'); hold on;
            xlabel('$\frac{|\Omega|}{2d(r-1)}$','Interpreter','latex')
            ylabel('evaluation')           
            subplot(1,2,2);
            plot(probiter*(d1*d2)/((d1+d2)*r*log(d1+d2)),...
                    mean(squeeze(resultbest(:,gi,i,:,k))),'-o'); hold on;
            xlabel('$\frac{|\Omega|}{2drlog(d)}$','Interpreter','latex')
            ylabel('evaluation')           

        end
    end
    legend({'Spearman Rho','Kendall Tau','RMSE'})
end

for gi=1:2:length(giter)
    figure()
    for i=1:length(dim_iter)
        [d1,d2,r]=dim_iter{i}{:};        
        for k=1:length(f)
            if k~=3
                temp=squeeze(resultbest(:,gi,i,:,k)>=0.999);               
            else
                temp=squeeze(resultbest(:,gi,i,:,k)<=1e-3);               
            end
            subplot(1,2,1);
            plot(probiter*(d1*d2)/((d1+d2-r)*r),...
                    mean(temp),'-o'); hold on;
            xlabel('$\frac{|\Omega|}{(2d-r)r}$','Interpreter','latex')
            ylabel('fraction of success')                     
            subplot(1,2,2);
            plot(probiter*(d1*d2)/((d1+d2)*r*log(d1+d2)),...
                    mean(temp),'-o'); hold on;
            xlabel('$\frac{|\Omega|}{2drlog(d)}$','Interpreter','latex')
            ylabel('fraction of success')                      
        end
    end
    title('Success=>ktau>0.999,srho>0.999,rmse<0.001, averaged over 30 runs')  
    legend({'Spearman Rho','Kendall Tau','RMSE'})
end
