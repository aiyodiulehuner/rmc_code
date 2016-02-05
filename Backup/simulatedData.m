%Simulating Data for RMC
function data=simulatedData(d,r)
if length(d)==1
    d1=d;d2=d;
else
    d1=d(1);d2=d(2);
end

data.U=randn(d1,r);%U=normc(U);
data.V=randn(d2,r);%V=normc(V);
data.rk=r;
data.Omega=rand(d1,d2);

theta=data.U*data.V';
data.mu=sum(svd(theta));
data.emin=inf;
for j=1:d2
    y=min(diff(sort(theta(:,j),'ascend')));
    if (y<data.emin)
        data.emin=y;
    end
end

    
%for b=0.1:0.1:1
%    a=[1:1:d]';
%    a=1./a;
%    a=a.^b;
%    D = diag(a);
%    U=randn(d,r);
%    V=randn(d,r);
%    U=D*U;
%    V=D*V;
%    data.U=U;
%    data.V=V;
%    data.rk=r;
%    data.Omega=rand(d,d);
%    disp([num2str(d) ' ' num2str(b)])
%    save(sprintf('simulatedDataUV_%d_%d.mat',d,round(10*b)),'data')   
