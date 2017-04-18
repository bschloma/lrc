
% Program:  lesEndOnly.m
% 
% Summary:  Compute sample paths over time of LES model but only save end points.  
%           Useful for computing stationary distributions.  
%
%           Populations grow logistically but are subject to multiplicatively coupled brownian motion.  
% 
% Inputs
%           r      : Bacterial growth rate (default 1.5 hr^{-1})
%           Kparams : Bacterial carrying capacity .
%                     Can be a single number, or [mean K, sigma of *log10*K for a
%                     log-normal distribution]. In the latter case, each run will
%                     draw K from a log-normal distribution. Default [1e4 0.5]
%           sig       : Volatility (default .8 hr^{-1/2}
%           numtrials  : number of sirlation runs
%           Tmax    : hrs., end time of sirlation (default 24 hours)
%           dt      : hrs., time step of sirlation (default 0.1 hours)
%           lextinct: Logical for extinction, true=popthresh (fixed at 1)
%                     is an absorbing state, false = no extinction
% 
% 
% Outputs
%           popend   : Array of final populations; columns = trials; 
%           
% 
% Fixed parameters
%           Initial population size : 10
% 
% Author:   Brandon Schlomann, modified from stochastic_collapse_model by Raghuveer Parthasarathy
% 
% Date:     3/18/17

function [popend] = lesEndOnly(r,Kparams,sig, numtrials, Tmax, dt,lextinct)

rng('shuffle'); %initialize random number generator

% Default parameter values

if ~exist('r', 'var') || isempty(r)
    r = 1.5; % 1/hours
end
if ~exist('Kparams', 'var') || isempty(Kparams)
    Kparams = [1e4 0.5]; % Carrying capacity, and sigma of log-normal distribution
end
if ~exist('Tmax', 'var') || isempty(Tmax)
    Tmax = 24; % hours
end
if ~exist('dt', 'var') || isempty(dt)
    dt = 0.1; % hours
end

if ~exist('numtrials', 'var') || isempty(numtrials)
    numtrials = 100;
end

if ~exist('sig', 'var') || isempty(sig)
    sig = .8;
end
if length(Kparams) == 1
    K = Kparams*ones(1,numtrials); % carrying capacity
else
    meanlogK = log10(Kparams(1));
    log10K = meanlogK + Kparams(2)*randn(1,numtrials);
    K = 10.^log10K;
end

if ~exist('lextinct', 'var') || isempty(lextinct)
    lextinct = false;
end



% Set things up
tarray = 0:dt:Tmax;  % array of time points
numsteps = length(tarray); % number of time points, including t0


%Initialize arrays
%popmat = zeros(numsteps,numtrials);
%popmat(1,:) = 10; % initial population size

popthresh = 1; 
pop0 = 10;

pop_prior = pop0*ones(1,numtrials);

% The actual simulation
for j=2:numsteps
     
    random_numbers = randn(1,numtrials); % pre-calculate random numbers from standard normal
    
    % Milstein method
    thispop = pop_prior + dt*r*pop_prior.*(1-pop_prior./K) + sig.*pop_prior.*random_numbers.*sqrt(dt) + .5.*pop_prior.*sig.^2.*dt.*(random_numbers.^2-1);
    
    if lextinct
        thispop(thispop<popthresh) = 0;
    end
    
    
    pop_prior = thispop;
    
    
end

popend = thispop;

end




