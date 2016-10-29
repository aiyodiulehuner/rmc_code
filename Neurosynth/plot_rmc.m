iters=[1,2,3,4,5];
probiter=0.2:0.2:0.8;

disp('ndgc,rho,tau')

load(sprintf('../results/ns_resultSMC.mat'))
val_results=squeeze(resultSMC(:,:,:,2,:));
val_results_cvmean=squeeze(mean(val_results,1));
msSMC=zeros(length(probiter),4,2);
niter=[1,2,3,4,5]
for pi=1:length(probiter)
    fprintf('p=%0.2g\n',probiter(pi))
    [~,ix_ndcg]=max(val_results_cvmean(pi,:,3));
    [~,ix_tau]=max(val_results_cvmean(pi,:,2));
    [~,ix_rho]=max(val_results_cvmean(pi,:,1));

     fprintf('\t ndcg=%d,rho=%d,tau=%d\n',ix_ndcg,ix_rho,ix_tau)
 
     msSMC(pi,1,1)=mean(resultSMC(:,pi,ix_ndcg,3,3));
     msSMC(pi,2,1)=mean(resultSMC(:,pi,ix_rho,3,1));
     msSMC(pi,3,1)=mean(resultSMC(:,pi,ix_tau,3,2));

     msSMC(pi,1,2)=std(resultSMC(:,pi,ix_ndcg,3,3));
     msSMC(pi,2,2)=std(resultSMC(:,pi,ix_rho,3,1));
     msSMC(pi,3,2)=std(resultSMC(:,pi,ix_tau,3,2)); 
     
     tp=zeros(length(niter),1);
     for i=1:length(niter)
         load(sprintf('../results/ns_ysmc_p%d/ysmc_cv%d_mi%d.mat',round(probiter(pi)*100),niter(i),ix_ndcg))
         tp(i)=t;
     end
     msSMC(pi,4,1)=mean(tp);msSMC(pi,4,2)=std(tp);
end

%clear resultSMC,val_results,val_results_cvmean,pi,ix_ndcg,ix_tau,ix_rho;

load(sprintf('../results/ns_resultRMC.mat'))
val_results=squeeze(resultRMC(:,:,:,2,:));
val_results_cvmean=squeeze(mean(val_results,1));
msRMC=zeros(length(probiter),4,2);

for pi=1:length(probiter)
    fprintf('p=%0.2g\n',probiter(pi))

    [~,ix_ndcg]=max(val_results_cvmean(pi,:,3));                                
    [~,ix_tau]=max(val_results_cvmean(pi,:,2));                                  
    [~,ix_rho]=max(val_results_cvmean(pi,:,1));                                  

     fprintf('\t ndcg=%d,rho=%d,tau=%d\n',ix_ndcg,ix_rho,ix_tau) 
     msRMC(pi,1,1)=mean(resultRMC(:,pi,ix_ndcg,3,3));                                          
     msRMC(pi,2,1)=mean(resultRMC(:,pi,ix_rho,3,1));                                          
     msRMC(pi,3,1)=mean(resultRMC(:,pi,ix_tau,3,2));                                          
                                                                                              
     msRMC(pi,1,2)=std(resultRMC(:,pi,ix_ndcg,3,3));                                           
     msRMC(pi,2,2)=std(resultRMC(:,pi,ix_rho,3,1));                                           
     msRMC(pi,3,2)=std(resultRMC(:,pi,ix_tau,3,2));               


     tp=zeros(length(niter),1);
     for i=1:length(niter)
         load(sprintf('../results/ns_yrmc_p%d/ns_yrmc_cv%d_mi%d.mat',round(probiter(pi)*100),niter(i),ix_ndcg))
         tp(i)=t;
     end
     msRMC(pi,4,1)=mean(tp);msRMC(pi,4,2)=std(tp);
     %load(sprintf('rmc%d/resultRMC.mat',i))
    %val_results=squeeze(resultRMC(1,:,:,3,:));
    %max from NDCG
end


close all
pli={'NDCG@10','Spearman Rho','Kendall Tau'};
np=length(pli);

csvwrite('../results/msRMC.csv',msSMC)
csvwrite('../results/msSMC.csv',msRMC)

% FIGURES IN plot_rmc.ipynb

h=figure();
h.set('Resize','off');
for i=1:np
disp(i)
subplot(1,np,i);
errorbar(probiter,msSMC(:,i,1),msSMC(:,i,2),'-bo','Linewidth',2);hold on
errorbar(probiter,msRMC(:,i,1),msRMC(:,i,2),'-ro','Linewidth',2)
xlim([0,1])
ylabel(pli(i));
xlabel('Fraction of matrix used in training')
end
set(h,'Position',[0,0,0.5,0.5])
