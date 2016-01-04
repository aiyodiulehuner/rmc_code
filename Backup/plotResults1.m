close all;clear all;clc
lamda=[1e-6,1e-5,1e-4,0.001,0.01,0.1,1,10,100,1000,1e4,1e5]';
d=1000;r=10;
g='argsort';
params=3:12;
SMCparam=lamda(params);

ndc=0;
ppfigSMC=0;

bb=2:2:10;
tau=zeros(length(bb),length(SMCparam));
rho=zeros(length(bb),length(SMCparam));
if (ndc)
        ndcg(p)=ndcg(U*V',data.estU,data.estV);
end


for i=1:length(bb);
    b=bb(i);
    load(sprintf('simulatedDataUV_%d_%d.mat',d,b));
    U=data.U;V=data.V;
    clear data

    for p=1:length(params)    
        load(sprintf('resultsSMC_%s_%d_%d_%d.mat',g,d,b,params(p)));
        tau(i,p)=data.kktau;rho(i,p)=data.ssrho;
        if (ndc)
            ndcg(i,p)=ndcg(U*V',data.estU,data.estV);
        end
    end
end

yl=0;

if (ppfigSMC)
htau=figure();
axes('Parent',htau,'XMinorTick','on','XScale','log','FontSize',12, 'Position',[0.15 0.15 0.8 0.8]);
for i=1:length(bb)
    htau=semilogx(SMCparam,tau(i,:),'-o','LineWidth',2);
    xlabel('$\lambda\frac{|\Omega|}{d_1d_2}$','interpreter', 'latex','FontSize',14);
    ylabel('Kendall Tau','interpreter', 'latex','FontSize',14); hold on;
    ylim([yl,1]);
    xlim([1e-6,1e4]);
    legend({'b=0.2','b=0.4','b=0.6','b=0.8','b=1.0'},'interpreter','latex','FontSize',14)
end

hrho=figure();
axes('Parent',hrho,'XMinorTick','on','XScale','log','FontSize',12, 'Position',[0.15 0.15 0.8 0.8]);
for i=1:length(bb)
    hrho=semilogx(SMCparam,rho(i,:),'-o','LineWidth',2);
    xlabel('$\lambda\frac{|\Omega|}{d_1d_2}$','interpreter', 'latex','FontSize',14);
    ylabel('Spearman Rho','interpreter', 'latex','FontSize',14); hold on;
    ylim([yl,1]);
    xlim([1e-6,1e4]);
    legend({'b=0.2','b=0.4','b=0.6','b=0.8','b=1.0'},'interpreter','latex','FontSize',14)
end

if (ndc)
    hndcg=figure();
    axes('Parent',hndcg,'XMinorTick','on','XScale','log','FontSize',12, 'Position',[0.15 0.15 0.8 0.8]);
    for i=1:length(bb)
        hndcg=semilogx(SMCparam,ndcg(i,:),'-o','LineWidth',2);
        xlabel('$\lambda\frac{|\Omega|}{d_1d_2}$','interpreter', 'latex','FontSize',14);
        ylabel('NDCG','interpreter', 'latex','FontSize',14); hold on;
        ylim([yl,1]);
        xlim([1e-6,1e4]);
        legend({'b=0.2','b=0.4','b=0.6','b=0.8','b=1.0'},'interpreter','latex','FontSize',14)
    end
end
end;









baselineData=load('results_baseline');

bestParamId=2;
lm=SMCparam(bestParamId);%*(d^2)/(2*d*r*log(2*d));
hb_tau=figure();
axes('Parent',hb_tau,'XMinorTick','on','XScale','log','FontSize',12, 'Position',[0.15 0.15 0.8 0.8]);
hb_tau=plot(bb/10.0,tau(:,bestParamId),'-o','LineWidth',2); hold on;
plot(bb/10.0,baselineData.data.kktau','-o','LineWidth',2);
xlabel('Coherence measure $b$','interpreter', 'latex','FontSize',14);
ylabel('Kendall Tau','interpreter', 'latex','FontSize',14); hold on;
ylim([yl,1]);
xlim([0,1]);
legend({['SMC ($\lambda=$' num2str(lm) ')'],'RMC ($\lambda=$)','Baseline'},'interpreter','latex','FontSize',14)

hb_rho=figure();
axes('Parent',hb_rho,'XMinorTick','on','XScale','log','FontSize',12, 'Position',[0.15 0.15 0.8 0.8]);
hb_rho=plot(bb/10.0,rho(:,bestParamId),'-o','LineWidth',2); hold on;
plot(bb/10.0,baselineData.data.ssrho','-o','LineWidth',2);
xlabel('Coherence measure $b$','interpreter', 'latex','FontSize',14);
ylabel('Spearman Rho','interpreter', 'latex','FontSize',14); hold on;
ylim([yl,1]);
xlim([0,1]);
legend({'SMC ($\lambda=$)','RMC ($\lambda=$)','Baseline'},'interpreter','latex','FontSize',14)