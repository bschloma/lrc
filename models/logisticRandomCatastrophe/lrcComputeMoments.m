% Program:  lrcComputeSamplePaths.m
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

function [moments,params] = lrcComputeMoments(runname,params)

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
lambda = params.lambda;              
lextinct = params.lextinct;       
M = params.M;               


%% Compute

[moments] = lrcMoments(mu,Kparams,f, lambda, numtrials, Tmax, dt,lextinct,M);
tvec = 0:dt:Tmax;

%% Plot
if lplot
    
    mvec = 1:M;
    scaled_moments = scaleByPowersOfK(moments,Kparams(1));
    
    figure; hold on;
    plot(mu.*tvec,scaled_moments,'linewidth',3);   
    set(gca,'fontsize',24,'linewidth',4)
    xlabel('rt','fontsize',24);
    ylabel('E[X^n_t/K^n]','fontsize',24);
    %axis([0 Tmax 0 .5);
    
    legendcell = cell(1,M);
    for mm =1:M
        legendcell{mm} = ['n = ' num2str(mm)];
    end
    
    legend(legendcell,'fontsize',18,'location','se');
  
end