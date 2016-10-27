clear;clc;clear globa;
for iX=1:5
outtest=sprintf('../../data/ml-100k/folds/cofirank_test_%d.lsvm',iX);
outval=sprintf('../../data/ml-100k/folds/cofirank_val_%d.lsvm',iX);
outtrain=sprintf('../../data/ml-100k/folds/cofirank_train_%d.lsvm',iX);
load(sprintf('../../data/ml-100k/folds/ml_%d.mat',iX));

file_train = fopen(outtrain,'w');
file_val = fopen(outval,'w');
file_test = fopen(outtest,'w');

for j=1:length(Jcol)-1
    disp(j)
    ind_train=Jcol(j)+1:Jcol(j+1);
    ind_val=Jcol_val(j)+1:Jcol_val(j+1);
    ind_test=Jcol_test(j)+1:Jcol_test(j+1);
    
    trainline='';
    if (~isempty(ind_train))
        ij=ii(ind_train);yj=yy(ind_train);
        [iis,idxs]=sort(ij,'ascend');yys=yj(idxs);    
        for i=1:length(iis)
            trainline=sprintf('%s%d:%d ',trainline,iis(i),yys(i));
        end
        trainline(end)=sprintf('\n');    
    else
        trainline='\n';
    end
    
    valline='';
    if (~isempty(ind_val))
        ij_val=ii_val(ind_val);yj_val=yy_val(ind_val);
        [ii_vals,idx_vals]=sort(ij_val,'ascend');yy_vals=yj_val(idx_vals);    
        for i=1:length(ii_vals)
            valline=sprintf('%s%d:%d ',valline,ii_vals(i),yy_vals(i));
        end
        valline(end)=sprintf('\n');
    else
        valline='\n';
    end
    
    testline='';
    if (~isempty(ind_test))
        ij_test=ii_test(ind_test);yj_test=yy_test(ind_test);
        [ii_tests,idx_tests]=sort(ij_test,'ascend');yy_tests=yj_test(idx_tests);       
        for i=1:length(ii_tests)
            testline=sprintf('%s%d:%d ',testline,ii_tests(i),yy_tests(i));
        end
        testline(end)=sprintf('\n');
    else
        testline='\n';
    end
    
    fprintf(file_train,trainline);
    fprintf(file_val,valline);
    fprintf(file_test,testline);
end
fclose(file_train);
fclose(file_val);
fclose(file_test);
end
