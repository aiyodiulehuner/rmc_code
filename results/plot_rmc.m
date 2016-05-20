iters=[1,2,3];
probiter=0.2:0.2:0.8;
testSMC=zeros(length(probiter),4,length(iters));
for ii=1:length(iters)
    i=iters(ii);
    load(sprintf('smc%d/resultSMC.mat',i))
    val_results=squeeze(resultSMC(1,:,:,3,:));
    %max from NDCG
    [~,ix_ndcg]=max(val_results(:,:,3),[],2);disp(ix_ndcg')
    [~,ix_tau]=max(val_results(:,:,1),[],2);disp(ix_tau')
    [~,ix_prec]=max(val_results(:,:,5),[],2);disp(ix_ndcg')


    for pi=1:length(probiter)
        testSMC(pi,1,ii)=resultSMC(1,pi,ix_ndcg(pi),3,1);
        testSMC(pi,2,ii)=resultSMC(1,pi,ix_tau(pi),3,3);
        testSMC(pi,3,ii)=resultSMC(1,pi,ix_prec(pi),3,5);
        load(sprintf('smc%d/ysmc_%d_%d.mat',i,round(probiter(pi)*100),ix_ndcg(pi)))
        testSMC(pi,4,ii)=t;
    end
end
mSMC=mean(testSMC,3);
eSMC=squeeze(std(permute(testSMC,[3,1,2])));


iters=[1];
probiter=0.2:0.2:0.8;
testRMC=zeros(length(probiter),4,length(iters));
for ii=1:length(iters)
    i=iters(ii);
    load(sprintf('rmc%d/resultRMC.mat',i))
    val_results=squeeze(resultRMC(1,:,:,3,:));
    %max from NDCG
    [~,ix_ndcg]=max(val_results(:,:,3),[],2);disp(ix_ndcg')
    [~,ix_tau]=max(val_results(:,:,1),[],2);disp(ix_tau')
    [~,ix_prec]=max(val_results(:,:,5),[],2);disp(ix_ndcg')


    for pi=1:length(probiter)
        testRMC(pi,1,ii)=resultRMC(1,pi,ix_ndcg(pi),3,1);
        testRMC(pi,2,ii)=resultRMC(1,pi,ix_tau(pi),3,3);
        testRMC(pi,3,ii)=resultRMC(1,pi,ix_prec(pi),3,5);
        load(sprintf('rmc%d/yrmc_%d_%d.mat',i,round(probiter(pi)*100),ix_ndcg(pi)))
        testRMC(pi,4,ii)=t;
    end
end
mRMC=mean(testRMC,3);
%eRMC=squeeze(std(permute(testRMC,[3,1,2])));
eRMC=zeros(size(mRMC));

probiter=0.2:0.2:0.8;
t=[10000,10000,10000,10000];
mCOFI=zeros(length(probiter),4);
load(sprintf('cofirank%d/resultCOFI.mat',i))
mCOFI(:,1:3)=squeeze(resultCOFI(1,1,:,3,[1,3,5]));
mCOFI(:,4)=t;
eCOFI=zeros(size(mCOFI));

close all
pli=[1,2,3];
np=length(pli);
for i=1:np
subplot(1,np,i); 
errorbar(probiter,mSMC(:,i),eSMC(:,i),'-o','Linewidth',2);hold on
errorbar(probiter(1:end),mRMC(:,i),eRMC(:,i),'-o','Linewidth',2)
xlim([0,1])
%ylim([0.5,1])
end

clear
niter=2;
testSMC=zeros(4,niter);
testRMC=zeros(4,niter);
testCOFI=zeros(4,niter);

for i=1:niter
    load(sprintf('smcml%d/ml_resultSMC.mat',i))
    val_results=squeeze(resultSMC(1,:,:,2,:));
    %max from NDCG
    [~,ix_ndcg]=max(val_results(:,3),[],1);disp(ix_ndcg')
    [~,ix_tau]=max(val_results(:,1),[],1);disp(ix_tau')
    [~,ix_prec]=max(val_results(:,5),[],1);disp(ix_ndcg')
    testSMC(1,i)=resultSMC(1,1,ix_ndcg,3,1);
    testSMC(2,i)=resultSMC(1,1,ix_tau,3,3);
    testSMC(3,i)=resultSMC(1,1,ix_prec,3,5);
    
    load(sprintf('rmcml%d/ml_resultRMC.mat',i))
    val_results=squeeze(resultRMC(1,:,:,2,:));
    %max from NDCG
    [~,ix_ndcg]=max(val_results(:,3),[],1);disp(ix_ndcg')
    [~,ix_tau]=max(val_results(:,1),[],1);disp(ix_tau')
    [~,ix_prec]=max(val_results(:,5),[],1);disp(ix_ndcg')
    testRMC(1,i)=resultRMC(1,1,ix_ndcg,3,1);
    testRMC(2,i)=resultRMC(1,1,ix_tau,3,3);
    testRMC(3,i)=resultRMC(1,1,ix_prec,3,5);
end
mSMC=mean(testSMC,2);
eSMC=squeeze(std(permute(testSMC,[2,1])));

mRMC=mean(testRMC,2);
eRMC=squeeze(std(permute(testRMC,[2,1])));

load('cofirank_ml50/ml_resultCOFI.mat')
testCOFI(1,1)=resultCOFI(1,1,1,3,1);
testCOFI(2,1)=resultCOFI(1,1,1,3,3);
testCOFI(3,1)=resultCOFI(1,1,1,3,5);

load('cofirank_ml50b/ml_resultCOFI.mat')
testCOFI(1,2)=resultCOFI(1,1,1,3,1);
testCOFI(2,2)=resultCOFI(1,1,1,3,3);
testCOFI(3,2)=resultCOFI(1,1,1,3,5);

mCOFI=mean(testCOFI,2);
eCOFI=squeeze(std(permute(testCOFI,[2,1])));
