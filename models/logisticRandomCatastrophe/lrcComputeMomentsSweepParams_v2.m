% Program:  lrcComputeMomentsSweepParams.m
%
% Summary:  Function to organize computation of LRC moments, sweeping both lambda
%           and f such that lambda*log(f) is constant for each pair of (lambda,f).
%
% Usage:    [momentsL,momentsF,params] = lrcComputeMomentsSweepParams(runname,params);
% 
% Inputs:   runname - string, name of computation for identification.
%                     Default = date and time string
%
%           params - struct, contains parameter values.  Default  = see
%                    below.  Can pass previously computed struct.
%
% Outputs:  momentsL - array containing moments from sweeping lambda.  
%                       rows = time, columns = moment #, depth = different lambda values
%           momentsF - array containing moments from sweeping f.
%                       rows = time, columns = moment #, depth = different f values
%           params -    struct, contains parameter values used in calculation
%
% Author:   Brandon Schlomann
%
% Date:     4/11/17 - first written, based off of older functions,
%                     including makeMeanVarZSweepPlot.m
%

function [momentsL,momentsF,params] = lrcComputeMomentsSweepParams_v2(runname,params)

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
   
    % struct to keep tract of params
    params = lrcParamsClass(runname);
    
    params = params.initZSweep();
%     
%     params.mu = 1;                 
%     params.meanK = 10^(4);         
%     params.sigK = 0;               
%     params.Kparams = [params.meanK params.sigK]; 
%     params.Tmax = 100;             
%     params.dt = .1;                
%     params.numtrials = 5e2;        
%     params.zmin = 1e-3;            
%     params.zmax = .8;              
%     params.numzs = 6;             
%     params.f = 1e-2;             
%     params.l = 1e-1;              
%     params.lextinct = false;       
%     params.M = 8;                  
%     params.zmin_an = 1e-3;
%     params.zmax_an = 1.1;
%     params.numzs_an = 300;
%     params.numsteps = length([0:params.dt:params.Tmax]);
%     params.lplot = true;
%     params.poiscolor =[.4,.1,.8];
%     params.bcolor = [.1,.7,.2];
    
end

%% Unpack params locally
mu = params.mu;                 
%meanK = params.meanK;        
%sigK = params.sigK;               
Kparams = params.Kparams; 
Tmax = params.Tmax;             
dt = params.dt;                
numtrials = params.numtrials;
lplot = params.lplot;
poiscolor =params.poiscolor;
bcolor = params.bcolor ;
f = params.f;              
lambda = params.l;              
lextinct = params.lextinct;       
M = params.M;               

zmin = params.zsweep.zmin;            
zmax = params.zsweep.zmax;              
numzs = params.zsweep.numzs;             
   
zmin_an = params.zsweep.zmin_an;
zmax_an = params.zsweep.zmax_an;
numzs_an = params.zsweep.numzs_an;
numsteps = length(0:params.dt:params.Tmax);

%% Arrays
zarray = linspace(zmin,zmax,numzs);
larray = -zarray./log(f);
farray = exp(-zarray./lambda);
zarray_an = linspace(zmin_an,zmax_an,numzs_an);
larray_an = -zarray_an./log(f);
farray_an = exp(-zarray_an./lambda);
momentsL = zeros(numsteps,M,numel(larray));
momentsF = zeros(numsteps,M,numel(larray));
statMomL = zeros(numzs_an,M);
statMomF = zeros(numzs_an,M);

%% Compute
% sweep lambda
for j = 1:numzs
    disp(['j = ' num2str(j)])
    
    momentsL(:,:,j) = lrcMoments(mu,Kparams,f,larray(j),numtrials,Tmax,dt,lextinct,M);
    
end
    
% sweep f      
for k = 1:numzs
    disp(['k = ' num2str(k)])

    momentsF(:,:,k) = lrcMoments(mu,Kparams,farray(k),lambda,numtrials,Tmax,dt,lextinct,M);

end


%% Plot
if lplot
    
    meanEndL = reshape(momentsL(end,1,:),1,numzs);
    varEndL = reshape(momentsL(end,2,:) - momentsL(end,1,:).^2,1,numzs);
    
    meanEndF = reshape(momentsF(end,1,:),1,numzs);
    varEndF = reshape(momentsF(end,2,:) - momentsF(end,1,:).^2,1,numzs);
    
    for jj = 1:numzs_an
        statMomL(jj,:) = lrcExactStationaryMoments(mu,Kparams(1),larray_an(jj),f,M);
    end
    
    for kk = 1:numzs_an
        statMomF(kk,:) = lrcExactStationaryMoments(mu,Kparams(1),lambda,farray_an(kk),M);
    end
    
    meanEndL_an = statMomL(:,1);
    meanEndF_an = statMomF(:,1);
    
    varEndL_an = statMomL(:,2) - meanEndL_an.^2;
    varEndF_an = statMomF(:,2) - meanEndF_an.^2;
   

    figure; hold on;
    plot(meanEndL./Kparams(1),varEndL./Kparams(1)./Kparams(1),'kd','markersize',24,'markerfacecolor',poiscolor);
    plot(meanEndF./Kparams(1),varEndF./Kparams(1)./Kparams(1),'ko','markersize',24,'markerfacecolor',poiscolor);
    plot(meanEndL_an./Kparams(1),varEndL_an./Kparams(1)./Kparams(1),'k-','linewidth',4)
    plot(meanEndF_an./Kparams(1),varEndF_an./Kparams(1)./Kparams(1),'k-','linewidth',4)
    
    set(gca,'fontsize',24,'linewidth',4)
    xlabel('E[X]/\kappa','fontsize',24);
    ylabel('Var[X]/\kappa^2','fontsize',24);
    axis([0 1 0 .3])

end