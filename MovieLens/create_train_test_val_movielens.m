clear;clc;clear global
addpath('../utils/');
d2=943;d1=1682;
warning('off')
dropped_user=[];
for i=1:5
    X=importdata(sprintf('../../data/ml-100k/folds/u%d.base',i));
    for u=1:d2
        if length(find(X(:,1)==u))<10
        	dropped_user(length(dropped_user)+1)=u;
		end
	end
end

uid=[];
Xtrain={};Xtest={};
for i=1:5
	X=importdata(sprintf('../../data/ml-100k/folds/u%d.base',i));
    X=X(~ismember(X(:,1),dropped_user),:);
    Xt=importdata(sprintf('../../data/ml-100k/folds/u%d.test',i));
    Xt=Xt(~ismember(Xt(:,1),dropped_user),:);
	Xtrain{i}=X;Xtest{i}=Xt;
    uid=unique([uid;X(:,1)]);
end
clc;
disp(dropped_user)
d2=length(uid);
rev_uid(uid)=1:length(uid);
    
for i=1:5
    disp(i)
    X=Xtrain{i};
    Xt=Xtest{i};

	X(:,1)=rev_uid(X(:,1));
    Xt(:,1)=rev_uid(Xt(:,1));

    disp([size(X),size(Xt)]);

	omega=rand(length(X),1);
	
	jj=X(omega<0.9,1);
	ii=X(omega<0.9,2);
	yy=X(omega<0.9,3);    
	[yy,ii,Jcol]=processInput(ii,jj,yy,d2);
	fprintf('\tTrain:%d, %d, %d\n', length(Jcol)-1,length(yy),length(unique(ii)))


	jj_val=X(omega>=0.9,1);
	ii_val=X(omega>=0.9,2);
	yy_val=X(omega>=0.9,3);    
	[yy_val,ii_val,Jcol_val]=processInput(ii_val,jj_val,yy_val,d2);
	fprintf('\tValidation:%d, %d, %d\n', length(Jcol_val)-1,length(yy_val),length(unique(ii_val)))

	jj_test=Xt(:,1);
	ii_test=Xt(:,2);
	yy_test=Xt(:,3);    
	[yy_test,ii_test,Jcol_test]=processInput(ii_test,jj_test,yy_test,d2);
	fprintf('\tTest: %d, %d, %d\n', length(Jcol_test)-1,length(yy_test),length(unique(ii_test)))

	save(sprintf('../../data/ml-100k/folds/ml_%i.mat',i),...
        'yy','yy_val','yy_test','ii','ii_val','ii_test',...
        'Jcol','Jcol_val','Jcol_test','jj','jj_val','jj_test','rev_uid','d1','d2')

end
