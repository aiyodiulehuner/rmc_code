clear;clc;clear global
addpath('utils/');
d2=943;d1=1682;

X=importdata(sprintf('../ml-100k/folds/ub.base'));
omega=rand(length(X),1);

jj=X(omega<0.9,1);
ii=X(omega<0.9,2);
yy=X(omega<0.9,3);    
[yy,ii,Jcol]=processInput(ii,jj,yy);
fprintf('\tTrain:%d, %d, %d\n', length(Jcol)-1,length(yy),length(unique(ii)))



jj_val=X(omega>=0.9,1);
ii_val=X(omega>=0.9,2);
yy_val=X(omega>=0.9,3);    
[yy_val,ii_val,Jcol_val]=processInput(ii_val,jj_val,yy_val);
fprintf('\tValidation:%d, %d, %d\n', length(Jcol_val)-1,length(yy_val),length(unique(ii_val)))

Xtest=importdata(sprintf('../ml-100k/folds/ub.test'));
jj_test=Xtest(:,1);
ii_test=Xtest(:,2);
yy_test=Xtest(:,3);    
[yy_test,ii_test,Jcol_test]=processInput(ii_test,jj_test,yy_test);
fprintf('\tTest: %d, %d, %d\n', length(Jcol_test)-1,length(yy_test),length(unique(ii_test)))

save(sprintf('../ml-100k/folds/ml_b.mat'),...
        'yy','yy_val','yy_test','ii','ii_val','ii_test',...
        'Jcol','Jcol_val','Jcol_test','jj','jj_val','jj_test','d1','d2')

