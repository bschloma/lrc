% computeExamplePaths.m
function [poispath,brownpath] = computeExamplePaths()

f = 1e-2; 
lambda = .7e-1; 
mu = 1;
K = 1e4;
Tmax = 48;
dt = .01;
relstd = 1e-7;
plotopt = false;
numtrials = 1;
p = 1;
lplot = true;

h = sqrt(-2*lambda*log(f));
%h = 1;

tvec = 0:dt:Tmax;

[poispath] = lrc(mu, K,f, lambda, numtrials, Tmax, dt, false);
[brownpath] = les(mu,K,h,numtrials,Tmax,dt,false);

if lplot
    pcolor =[.4,.1,.8];
    bcolor = [.1,.7,.2];
    figure; hold on;
    plot(tvec,poispath./K,'linewidth',5,'color',pcolor)
    plot(tvec,brownpath./K,'linewidth',5,'color',bcolor)
    %xlabel('Time (hrs)','fontsize',24)
    %ylabel('Population','fontsize',24)
    set(gca,'fontsize',36)
    set(gca,'yscale','log')
    axis([0 50 10^-4 3])
    set(gca,'ytick',[10^-2 10^0 ])
end

end

