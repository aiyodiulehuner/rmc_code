%Simulated Data
%||Theta||=\sqrt{d1d2}

clear;clc;%close all;

niter=5;
d1=1000;d2=1000;r=20;
c=1:5;

sI=zeros(niter,d1);
sexp=zeros(niter,d1);
slog=zeros(niter,d1);
slogit=zeros(niter,d1,length(c));
sargsort=zeros(niter,d1);
srandom=zeros(niter,d1);
glist={@(a) a, @exp, @(a) log(a-min(a)+1), @(a) rand(1)*a + rand(1), @(a,c) 1./(1+exp(-c*a)) };


for i=1:niter
    
    U=randn(d1,r);
    V=randn(d2,r);
    
    theta=U*V'; 
    theta=sqrt(d1*d2)*theta/norm(theta,'fro');
     
    s=svd(theta); sI(i,:)=sqrt(cumsum(s.^2))/norm(s,2); 
    
    %exp
    Y=exp(theta); s=svd(Y); sexp(i,:)=sqrt(cumsum(s.^2))/norm(s,2); 
    
    %log
    Y=theta;m=min(Y);
    for j=1:d2
        Y(:,j)=log(theta(:,j)-m(j)+1); 
    end
    s=svd(Y); slog(i,:)=sqrt(cumsum(s.^2))/norm(s,2); 
    
    %argsort
    Y=theta;
    for j=1:d2
        [y,IX]=sort(Y(:,j),'ascend');
        Y(IX,j)=1:length(IX);
    end
    s=svd(Y); sargsort(i,:)=sqrt(cumsum(s.^2))/norm(s,2); 
    
    %random
    Y=theta;
    for j=1:d2
        Y(:,j)=glist{randsample(4,1)}(Y(:,j));
    end
    s=svd(Y); srandom(i,:)=sqrt(cumsum(s.^2))/norm(s,2);     
    
    %logit    
    for ci=1:length(c)
        Y=theta;
        for j=1:d2
            Y(:,j)=glist{5}(Y(:,j),c(ci));
        end
        s=svd(Y); slogit(i,:,ci)=sqrt(cumsum(s.^2))/norm(s,2); 
    end
  
end

% figure()
% % errorbar(1:d1,mean(sI),std(sI),'-o');hold on;
% % errorbar(1:d1,mean(sexp),std(sexp),'-o');hold on;
% % errorbar(1:d1,mean(slog),std(slog),'-o');hold on;
% % errorbar(1:d1,mean(sargsort),std(sargsort),'-o'); hold on;
% % errorbar(1:d1,mean(srandom),std(srandom),'-o'); hold on;
% % errorbar(1:d1,mean(slogit),std(slogit),'-o'); hold on;
% semilogx(1:d1,mean(sI),'-o');hold on;
% semilogx(1:d1,mean(sexp),'-o');hold on;
% semilogx(1:d1,mean(slog),'-o');hold on;
% semilogx(1:d1,mean(sargsort),'-o');hold on;
% semilogx(1:d1,mean(srandom),'-o'); hold on;
% txt={'I','exp','log','argsort','random'};
% for ci=2:2:length(c)
%     semilogx(1:d1,mean(slogit(:,:,ci)),'-o'); hold on;
%     txt{length(txt)+1}=['logit_c=', num2str(c(ci))];
% end
% legend(txt)
% 
% figure()
reps=zeros(1,length(c));
for ci=1:length(c)
    reps(ci)=find(mean(slogit(:,:,ci))>0.99,1);
end
plot(c',reps','-o'); hold on 


% for q=1:length(qiter)
%     semilogx(1:1000,mean(squeeze(squant(:,q,:))),'-o','Linewidth',2); hold on;
%     txt{q+5}=sprintf('QuantLevel %d',qiter(q));
% end
% legend(txt)
% ylim([0,1])
% quant
%     for q=1:length(qiter)
%         Y=theta;
%         m=min(min(Y)); M=max(max(Y)); 
%         if mod(qiter(q),2)
%             for j=1:qiter(q)
%                 Y(Y>=m+((j-1)/qiter(q))*(M-m) & Y<m+(j/qiter(q))*(M-m))=-floor(qiter(q)/2)+j-1;
%             end
%             s=svd(Y); squant(i,q,:)=s/max(s); 
%         else
%             for j=1:qiter(q)/2
%                 Y(Y>=m+((j-1)/qiter(q))*(M-m) & Y<m+(j/qiter(q))*(M-m))=-(qiter(q)/2)+j-1;
%             end
%             for j=(qiter(q)/2+1):qiter(q)
%                 Y(Y>=m+((j-1)/qiter(q))*(M-m) & Y<m+(j/qiter(q))*(M-m))=-(qiter(q)/2)+j;
%             end
%         end   
%         s=svd(Y); squant(i,q,:)=s/max(s);
%     end    