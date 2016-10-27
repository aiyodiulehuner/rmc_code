clear;clc;
niter=5;
msSMC=zeros(4,2);
msRMC=zeros(4,2);
msCOFI=zeros(4,2);
disp(['ndcg','prec','tau','rho'])

%SMC
load(sprintf('../results/ml_resultSMC.mat'))
val_results=squeeze(resultSMC(:,:,2,:));
val_results_cvmean=squeeze(mean(val_results,1));

[~,ix_ndcg]=max(val_results_cvmean(:,3));disp(ix_ndcg')
[~,ix_tau]=max(val_results_cvmean(:,2));disp(ix_tau')
[~,ix_prec]=max(val_results_cvmean(:,5));disp(ix_prec')
[~,ix_rho]=max(val_results_cvmean(:,1));disp(ix_rho')

msSMC(1,1)=mean(resultSMC(:,ix_ndcg,3,3));
msSMC(2,1)=mean(resultSMC(:,ix_prec,3,5));
msSMC(3,1)=mean(resultSMC(:,ix_tau,3,2));
msSMC(4,1)=mean(resultSMC(:,ix_tau,3,1));

msSMC(1,2)=std(resultSMC(:,ix_ndcg,3,3));
msSMC(2,2)=std(resultSMC(:,ix_prec,3,5));                                                    
msSMC(3,2)=std(resultSMC(:,ix_tau,3,2));                                                     
msSMC(4,2)=std(resultSMC(:,ix_tau,3,1));

disp(msSMC)

%%RMC
%load(sprintf('../results/ml_resultRMC.mat'))
%val_results=squeeze(resultRMC(:,:,2,:));
%val_results_cvmean=squeeze(mean(val_results,1));
    
%[~,ix_ndcg]=max(val_results_cvmean(:,3));disp(ix_ndcg')
%[~,ix_tau]=max(val_results_cvmean(:,2));disp(ix_tau')
%[~,ix_prec]=max(val_results_cvmean(:,5));disp(ix_prec')
%[~,ix_rho]=max(val_results_cvmean(:,1));disp(ix_rho')

%msRMC(1,1)=mean(resultRMC(:,ix_ndcg,3,3));                                                    
%msRMC(2,1)=mean(resultRMC(:,ix_prec,3,5));                                                    
%msRMC(3,1)=mean(resultRMC(:,ix_tau,3,2));                                                     
%msRMC(4,1)=mean(resultRMC(:,ix_tau,3,1));                                                     
                                                                                              
%msRMC(1,2)=std(resultRMC(:,ix_ndcg,3,3));
%msRMC(2,2)=std(resultRMC(:,ix_prec,3,5)); 
%msRMC(3,2)=std(resultRMC(:,ix_tau,3,2)); 
%msRMC(4,2)=std(resultRMC(:,ix_tau,3,1));

%disp(msRMC)
%}
%COFIRANK
load('../results/ml_resultCOFI.mat')
msCOFI(1,1)=mean(resultCOFI(:,1,3,3));
msCOFI(2,1)=mean(resultCOFI(:,1,3,5));
msCOFI(3,1)=mean(resultCOFI(:,1,3,2));
msCOFI(4,1)=mean(resultCOFI(:,1,3,1));

msCOFI(1,2)=std(resultCOFI(:,1,3,3));
msCOFI(2,2)=std(resultCOFI(:,1,3,5));
msCOFI(3,2)=std(resultCOFI(:,1,3,2));
msCOFI(4,2)=std(resultCOFI(:,1,3,1));

disp(msCOFI)
