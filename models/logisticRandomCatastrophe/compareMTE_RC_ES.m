function [mte_rc,mte_es] = compareMTE_RC_ES(r,K,dt,lambda,f,numtrials)

sigma = sqrt(-2.*lambda.*log(f));

numls = numel(lambda);
mte_rc = zeros(numls,1);
mte_es = zeros(numls,1);

% m_rc = zeros(numls,1);
% m_es = zeros(numls,1);
% v_rc = zeros(numls,1);
% v_es = zeros(numls,1);

pcolor =[.4,.1,.8];
bcolor = [.1,.7,.2];

Tmax = 100;

for n = 1:numls
    disp(num2str(n))
    mte_rc(n) = lrcMTE(r,K,lambda(n),f,dt,numtrials,false);
%     pend_rc = logisticGrowthPoissonCollapse_fast(f,lambda(n),numtrials,Tmax,dt,r,K);
% 
%     m_rc(n) = mean(pend_rc);
%     v_rc(n) = var(pend_rc);
    
    mte_es(n) = lesMTE(r,K,sigma(n),dt,numtrials);
    
%     pend_es = logisticGrowthBrownian_endOnly(r,K,sigma(n),dt,Tmax,numtrials,1);
% 
%     m_es(n) = mean(pend_es);
%     v_es(n) = var(pend_es);
    
end

figure; hold on;
plot(sigma.^2/2,mte_rc,'ko','markersize',24,'markerfacecolor',pcolor)
plot(sigma.^2/2,mte_es,'ko','markersize',24,'markerfacecolor',bcolor)
set(gca,'fontsize',24)
set(gca,'yscale','log')
title('MTE','fontsize',24)
legend({'LRC','LES'},'fontsize',24)

% figure; hold on;
% plot(sigma.^2/2,m_rc./K,'ko','markersize',24,'markerfacecolor',pcolor)
% plot(sigma.^2/2,m_es./K,'ko','markersize',24,'markerfacecolor',bcolor)
% set(gca,'fontsize',24)
% title('Mean','fontsize',24)
% legend({'LRC','LES'},'fontsize',18)
% 
% figure; hold on;
% plot(sigma.^2/2,v_rc./K./K,'ko','markersize',24,'markerfacecolor',pcolor)
% plot(sigma.^2/2,v_es./K./K,'ko','markersize',24,'markerfacecolor',bcolor)
% set(gca,'fontsize',24)
% title('Variance','fontsize',24)
% legend({'LRC','LES'},'fontsize',18)



end