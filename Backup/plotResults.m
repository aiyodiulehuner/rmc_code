close all
lamda=[1e-6,1e-5,1e-4,0.001,0.01,0.1,1,10,100,1000,1e4,1e5]';
gg={'I','scale_shift','argsort','log','exp'};
param_iter={3:12,4:12,3:12,1:12,1:12};

h=figure();
axes('Parent',h,'XMinorTick','on','XScale','log','FontSize',12, 'Position',[0.15 0.15 0.8 0.8]);
for i=[3];
g=gg{i};
pi=param_iter{i};
SMCparam=lamda(pi);
tau=zeros(size(SMCparam));
rho=zeros(size(SMCparam));
rmse=zeros(size(SMCparam));
for p=1:length(pi)    
    load(sprintf('resultsSMC_%s_%d_%d.mat',g,d1,pi(p)));
    tau(p)=data.ktau;rho(p)=data.srho;
end
semilogx(SMCparam,tau,'-o','LineWidth',2);
xlabel('$\lambda\frac{|\Omega|}{d_1d_2}$','interpreter', 'latex','FontSize',14);
ylabel('Kendall Tau','interpreter', 'latex','FontSize',14); hold on;
%ylim([0,1]);
xlim([1e-6,1e4]);
end
legend({'$Y_{I(.)}$','$Y_{a(.)+b}$','$Y_{{argsort(.)}}$', '$Y_{\log{(.+b)}}$','$Y_{\exp(.)}$'},...
    'interpreter','latex','FontSize',14)

%save results

% select param
% print table result
% plot time result
