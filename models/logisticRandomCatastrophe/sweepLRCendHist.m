

function [pophistmat,popbins,params] = sweepLRCendHist(runname,params)

%% Ensure correct input
if ~exist('runname', 'var') || isempty(runname)
    % given name based on date and time
    c= clock;
    runname = '';
    for cc = 1:5   %leave off seconds
        runname = [runname num2str(c(cc)) '_'];
    end
    
    runname = runname(1:(end-1));    %remove trailing _
            
end


if ~exist('params', 'var') || isempty(params)
   
    % struct to keep track of params
    params = struct();
    params.runname = runname;
    
    params.mu = 1;                 
    params.meanK = 10^(4);         
    params.sigK = 0;               
    params.Kparams = [params.meanK params.sigK]; 
    params.Tmax = 50;             
    params.dt = .001;                
    params.numtrials = 1e5;        
    %params.zmin = 1e-3;            
    %params.zmax = .8;              
    %params.numzs = 6;        
    params.f0 = 1e-2;             
    params.l0 = 1e-1;              
    params.lextinct = false;       
    %params.M = 8;                  
    %params.zmin_an = 1e-3;
    %params.zmax_an = 1.1;
    %params.numzs_an = 300;
    params.lplot = true;
    params.poiscolor =[.4,.1,.8];
    params.bcolor = [.1,.7,.2];
    
    %params.sigmin = .0047;            
    %params.sigmax = 1.2649;              
    %params.numsigs = 6;      
    %params.sigmin_an = .0047;
    %params.sigmax_an = 1.2649;
    %params.numsigs_an = 300;
    
    params.numbins = 20;

    params.numalphs = 5;
    params.sigma = .8;
    params.lmultmax = 75;

    params.ymax = .3;
end

%% Unpack params locally
mu = params.mu;                               
Kparams = params.Kparams;
K = Kparams(1);
Tmax = params.Tmax;             
dt = params.dt;                
numtrials = params.numtrials;                
f0 = params.f0;              
l0 = params.l0;              
lextinct = params.lextinct;       
lplot = params.lplot;
numbins = params.numbins;
numalphs = params.numalphs;
sigma = params.sigma;
lmultmax = params.lmultmax;
poiscolor = params.poiscolor;
bcolor = params.bcolor;
ymax = params.ymax;

%% Arrays
scalevec = linspace(1,lmultmax,numalphs);
lvec = l0.*scalevec;
fvec = exp(-sigma./sqrt(lvec));
rvec = mu -lvec.*log(fvec) - lvec.*log(fvec).*log(fvec)./2;
Kvec = K.*(1-lvec.*log(fvec)./mu - lvec.*log(fvec).*log(fvec)./2./mu);

mendvec = zeros(1,numalphs);
varvec = zeros(1,numalphs);

pophistmat = zeros(numalphs,numbins);

if lplot
    figure; hold on;
end

%% Loop through alphas
for s = 1:numalphs
    %popbins = logspace(-10,log10(max(Kvec)),numbins);
    [pophist,popbins,mendvec(s),varvec(s)] = computeLRCendHist(rvec(s),Kvec(s),lvec(s),fvec(s),dt,Tmax,numtrials,lextinct,numbins,false);
    if lplot
        ax_s = subplot(2,3,s);
        bar(ax_s,popbins,pophist./sum(pophist),'facecolor',poiscolor,'edgecolor',poiscolor)
        %xlabel('log10(X)','fontsize',22)
        %ylabel('P(log10(X))','fontsize',22)
        set(gca,'fontsize',22)
        
        if s > 1
            axis([0 5 0 ymax])
        else
            axis([0 5 0 .8])
        end
    end
    
    pophistmat(s,:) = pophist;
    
end

 s = s+1;
 rB = mu - sigma^2/2;
 KB = K*(1-sigma^2/2/mu);
 dtB = .01;
 
 [pophist,popbins] = makeESendHist(rB,KB,sigma,dtB,Tmax,numtrials,numbins,false);
 ax_s = subplot(2,3,s);
 bar(ax_s,popbins,pophist./sum(pophist),'facecolor',bcolor,'edgecolor',bcolor)
 %xlabel('log10(X)','fontsize',22)
 %ylabel('P(log10(X))','fontsize',22)
 set(gca,'fontsize',22)
 axis([0 5 0 ymax])

% mend_theory = Kvec.*(1+lvec.*log(fvec)./rvec);
% figure; hold on;
% %plot(scalevec,mendvec,'ko','markerfacecolor',pcolor,'linewidth',2)
% plot(scalevec,mend_theory,'color',pcolor,'linewidth',4)
% set(gca,'fontsize',22)
% xlabel('\alpha','fontsize',22)
% ylabel('E[X]/K','fontsize',22)
% 
% vend_theory = Kvec.^2.*(-log(fvec)-(1-fvec)).*lvec./rvec.*(1+lvec.*log(fvec)./rvec);
% vend_theory_lim = Kvec.^2.*(log(fvec).^2).*lvec./2./rvec.*(1+lvec.*log(fvec)./rvec);
% figure; hold on;
% %plot(scalevec,varvec,'ko','markerfacecolor',pcolor,'linewidth',2)
% plot(scalevec,vend_theory./vend_theory_lim,'color',pcolor,'linewidth',4)
% plot(scalevec,vend_theory_lim./vend_theory_lim,'color','r','linewidth',4)
% set(gca,'fontsize',22)
% xlabel('\alpha','fontsize',22)
% ylabel('Var[X]','fontsize',22)
% set(gca,'xscale','log')
% axis([1 1e4 .5 1])

end

