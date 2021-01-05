function [] = f_plotposterior(out)
%F_PLOTPOSTERIOR Summary of this function goes here
%   Detailed explanation goes here

col = lines;

SYS(1,:) = [out.priors.SI(1)*1e4 out.priors.SI(2)];
SYS(2,:) = [out.priors.p1(1)*1e3 out.priors.p1(2)];
SYS(3,:) = [out.priors.p2(1)*1e3 out.priors.p2(2)];

SYS(4,:) = [out.posterior.SI(1)*1e4 out.posterior.SI(2)];
SYS(5,:) = [out.posterior.p1(1)*1e3 out.posterior.p1(2)];
SYS(6,:) = [out.posterior.p2(1)*1e3 out.posterior.p2(2)];

SYS(:,2) = exp(sqrt(log(SYS(:,2)/100+1)));

INP = [out.priors.k;out.posterior.k];


figure('Units','Normalized','Position',[0.05 0.5 0.4 0.3]),

subplot(121); hold on; box on;
title('System Parameters');

for i=1:3
    errorbar(i-0.1,SYS(i,1),SYS(i,1)*(1-1/SYS(i,2)),SYS(i,1)*(SYS(i,2)-1),...    
        'o','Color',col(1,:),'MarkerfaceColor',col(1,:),'LineWidth',1.2);
end
for i=4:6
    errorbar(i-3+0.1,SYS(i,1),SYS(i,1)*(1-1/SYS(i,2)),SYS(i,1)*(SYS(i,2)-1),...    
        'o','Color',col(2,:),'MarkerfaceColor',col(2,:),'LineWidth',1.2);
end

subplot(122); hold on; box on;
title(['Input Parameters for ' out.options.GA_fun]);

if strcmp(out.options.GA_fun,'RaPL')
    INP(:,2) = exp(sqrt(log(INP(:,2)/100+1)));
    for i=1:6
        errorbar(i-0.1,INP(i,1),INP(i,1)*(1-1/INP(i,2)),INP(i,1)*(INP(i,2)-1),...    
            'o','Color',col(1,:),'MarkerfaceColor',col(1,:),'LineWidth',1.2);
    end
    for i=7:12
        errorbar(i-6+0.1,INP(i,1),INP(i,1)*(1-1/INP(i,2)),INP(i,1)*(INP(i,2)-1),...    
            'o','Color',col(2,:),'MarkerfaceColor',col(2,:),'LineWidth',1.2);
    end
else
    INP([1:4,6:9],2) = exp(sqrt(log(INP([1:4,6:9],2)/100+1)));
    [INP,s_o] = f_logistic_mapping(m_i,s_i,mod)
    for i=1:6
        errorbar(i-0.1,INP(i,1),INP(i,1)*(1-1/INP(i,2)),INP(i,1)*(INP(i,2)-1),...    
            'o','Color',col(1,:),'MarkerfaceColor',col(1,:),'LineWidth',1.2);
    end
    for i=7:12
        errorbar(i-6+0.1,INP(i,1),INP(i,1)*(1-1/INP(i,2)),INP(i,1)*(INP(i,2)-1),...    
            'o','Color',col(2,:),'MarkerfaceColor',col(2,:),'LineWidth',1.2);
    end




end

