% Program:  lesComputeMomentsSweepParams.m
%
% Summary:  Function to organize computation of LES moments, sweeping
%           sigma.  
%
% Usage:    [momentsS,params] = lesComputeMomentsSweepParams(runname,params);
%
%           Can map sigma onto lambda, f values of LRC model by first calling
%           mapSigma.m with params = mapSigma(params);
%
% Inputs:   runname - string, name of computation for identification.
%                     Default = date and time string
%
%           params - struct, contains parameter values.  Default  = see
%                    below.  Can pass previously computed struct.
%
% Outputs:  momentsS - array containing moments from sweeping lambda.  
%                       rows = time, columns = moment #, depth = different sigma values
%           
%           params -    struct, contains parameter values used in calculation
%
% Author:   Brandon Schlomann
%
% Date:     4/11/17 - first written, based off of older functions,
%                     including makeMeanVarZSweepPlot.m
%

function [momentsS,params] = lesComputeMomentsSweepParams(runname,params)

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
    params = struct();
    params.runname = runname;
    
    params.mu = 1;                 
    params.meanK = 10^(4);         
    params.sigK = 0;               
    params.Kparams = [params.meanK params.sigK]; 
    params.Tmax = 100;             
    params.dt = .1;                
    params.numtrials = 5e2;        
    params.zmin = 1e-3;            
    params.zmax = .8;              
    params.numzs = 6;        
    params.f0 = 1e-2;             
    params.l0 = 1e-1;              
    params.lextinct = false;       
    params.M = 8;                  
    params.zmin_an = 1e-3;
    params.zmax_an = 1.1;
    params.numzs_an = 300;
    %params.numsteps = length([0:params.dt:params.Tmax]);
    params.lplot = true;
    params.poiscolor =[.4,.1,.8];
    params.bcolor = [.1,.7,.2];
    
    params.sigmin = .0047;            
    params.sigmax = 1.2649;              
    params.numsigs = 6;      
    params.sigmin_an = .0047;
    params.sigmax_an = 1.2649;
    params.numsigs_an = 300;
end

%% Unpack params locally
mu = params.mu;                 
meanK = params.meanK;        
sigK = params.sigK;               
Kparams = params.Kparams; 
Tmax = params.Tmax;             
dt = params.dt;                
numtrials = params.numtrials;        
sigmin = params.sigmin;            
sigmax = params.sigmax;              
numsigs = params.numsigs;             
f0 = params.f0;              
l0 = params.l0;              
lextinct = params.lextinct;       
M = params.M;                  
sigmin_an = params.sigmin_an;
sigmax_an = params.sigmax_an;
numsigs_an = params.numsigs_an;
numsteps = length([0:params.dt:params.Tmax]);
lplot = params.lplot;
poiscolor =params.poiscolor;
bcolor = params.bcolor ;
zmin_an = params.zmin_an;
zmax_an = params.zmax_an;
numzs_an = params.numzs_an;
zmin = params.zmin;            
zmax = params.zmax;              
numzs = params.numzs;   
%% Arrays
sigarray = linspace(sigmin,sigmax,numsigs);

sigarray_an = linspace(sigmin_an,sigmax_an,numsigs_an);

momentsS = zeros(numsteps,M,numel(sigarray));
statMomS = zeros(numsigs_an,M);

%% Compute
% sweep lambda
for j = 1:numsigs
    disp(['Iteration ' num2str(j) ', sigma = ' num2str(sigarray(j))])
    
    momentsS(:,:,j) = lesMoments(mu,Kparams,sigarray(j),numtrials,Tmax,dt,lextinct,M);
    
end
    

%% Plot
if lplot
    
    meanEndS = reshape(momentsS(end,1,:),1,numsigs);
    varEndS = reshape(momentsS(end,2,:) - momentsS(end,1,:).^2,1,numsigs);
    

    
    for jj = 1:numsigs_an
        statMomS(jj,:) = lesExactStationaryMoments(mu,Kparams(1),sigarray_an(jj),M);
    end
    
   
    
    meanEndS_an = statMomS(:,1);
    
    varEndS_an = statMomS(:,2) - meanEndS_an.^2;

    figure; hold on;
    plot(meanEndS./Kparams(1),varEndS./Kparams(1)./Kparams(1),'ks','markersize',24,'markerfacecolor',bcolor);
    plot(meanEndS_an./Kparams(1),varEndS_an./Kparams(1)./Kparams(1),'k-','linewidth',4)
    
    set(gca,'fontsize',24,'linewidth',4)
    xlabel('E[X]/K','fontsize',24);
    ylabel('Var[X]/K^2','fontsize',24);
    axis([0 1 0 .3])

end