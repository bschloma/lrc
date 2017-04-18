function [pophist,popbins,mend,vend] = computeLRCendHist(r,K,lambda,f,dt,Tmax,numtrials,lextinct,popbins,lplot)

[popend] = lrcEndOnly(r,K,f,lambda,numtrials,Tmax,dt,lextinct);

mend = mean(popend);
vend = var(popend);

[pophist,popbins] = hist(log10(popend),popbins);

if lplot
    figure; hold on;
    bar(popbins,pophist./(sum(pophist)))
    xlabel('X','fontsize',22)
    ylabel('P(X)','fontsize',22)
    set(gca,'fontsize',22)
end



end