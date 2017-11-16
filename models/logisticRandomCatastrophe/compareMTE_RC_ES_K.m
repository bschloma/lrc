function [mte_rc_K,mte_es_K,Kvec_LRC,Kvec_LES] = compareMTE_RC_ES_K(r,K,dt,lambda,f,numtrials,maptype,method,x0)


switch maptype
    case 1
        % equate stationary means
        sigma = sqrt(-2.*lambda.*log(f));
        
        rvec_LRC = r.*ones(1,numel(lambda));
        Kvec_LRC = K.*ones(1,numel(lambda));
        rvec_LES = rvec_LRC;
        Kvec_LES = Kvec_LRC;
    case 2
        % diffusion limit
        sigma = -sqrt(lambda).*log(f);
        rvec_LRC = r - lambda.*log(f);
        Kvec_LRC = K.*(1-lambda./r.*log(f));
        rvec_LES = r + sigma.^2./2;
        Kvec_LES = K.*(1+sigma.^2./r./2);
end
        
if ~exist('method', 'var') || isempty(method)
    method = 'pdmp';
end
if ~exist('x0', 'var') || isempty(x0)
    x0 = 10;
end


numKs = numel(K);
mte_rc_K = zeros(numKs,1);
mte_es_K = zeros(numKs,1);

%m_rc = zeros(numls,1);
%m_es = zeros(numls,1);
%v_rc = zeros(numls,1);
%v_es = zeros(numls,1);

pcolor =[.4,.1,.8];
bcolor = [.1,.7,.2];

Tmax = 100;

for n = 1:numKs
    disp(num2str(n))
    %mte_rc_K(n) = lrcMTE(r,K(n),lambda,f,dt,numtrials,false);
    mte_rc_K(n) = lrcMTE(rvec_LRC,Kvec_LRC(n),lambda,f,dt,numtrials,false,method,x0);

    %pend_rc = logisticGrowthPoissonCollapse_fast(f,lambda(n),numtrials,Tmax,dt,r,K);

%     m_rc(n) = mean(pend_rc);
%     v_rc(n) = var(pend_rc);
    
    mte_es_K(n) = lesMTE(rvec_LES,Kvec_LES(n),sigma,dt,numtrials,x0);
    
%     pend_es = logisticGrowthBrownian_endOnly(r,K,sigma(n),dt,Tmax,numtrials,1);
% 
%     m_es(n) = mean(pend_es);
%     v_es(n) = var(pend_es);
    
end

% figure; hold on;
% loglog(K,mte_rc_K,'ko','markersize',24,'markerfacecolor',pcolor)
% loglog(Kvec_LES,mte_es_K,'ks','markersize',24,'markerfacecolor',bcolor)
% set(gca,'fontsize',24)
% set(gca,'yscale','log')
% set(gca,'xscale','log')
% title('MTE','fontsize',24)
% legend({'LRC','LES'},'fontsize',24,'location','nw')

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
