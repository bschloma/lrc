
% Program:  lesMoments.m
% 
% Summary:  Compute moments over time of LES model.  Populations grow logistically but are
%           subject to rltiplicatively coupled brownian motion.  
% 
% Inputs
%           r      : Bacterial growth rate (default 1.5 hr^{-1})
%           Kparams : Bacterial carrying capacity .
%                     Can be a single number, or [mean K, sigma of *log10*K for a
%                     log-normal distribution]. In the latter case, each run will
%                     draw K from a log-normal distribution. Default [1e4 0.5]
%           h       : Volatility (default .8 hr^{-1/2}
%           numtrials  : number of sirlation runs
%           Tmax    : hrs., end time of sirlation (default 24 hours)
%           dt      : hrs., time step of sirlation (default 0.1 hours)
%           lextinct: Logical for extinction, true=popthresh (fixed at 1)
%                     is an absorbing state, false = no extinction
% 
% 
% Outputs
%           moments   : Array of moments; row = time; columns = trials; 
%           
% 
% Fixed parameters
%           Initial population size : 10
% 
% Author:   Brandon Schlomann, modified from stochastic_collapse_model by Raghuveer Parthasarathy
% 
% Date:     3/18/17

function [moments] = lesMoments(r,Kparams,h, numtrials, Tmax, dt,lextinct,M)

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

if ~exist('h', 'var') || isempty(h)
    h = .8;
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

if ~exist('M', 'var') || isempty(M)
    M = 6;
end

% Set things up
tarray = 0:dt:Tmax;  % array of time points
numsteps = length(tarray); % number of time points, including t0


%Initialize arrays
%popmat = zeros(numsteps,numtrials);
%popmat(1,:) = 10; % initial population size

popthresh = 1; 
mvec = 1:M;
moments = zeros(numsteps,M);
pop0 = 10;

pop_prior = pop0*ones(1,numtrials);
moments(1,:) = pop0.^mvec;

% The actual simulation
for j=2:numsteps
     
    random_numbers = randn(1,numtrials); % pre-calculate random numbers from standard normal
    
    % Milstein method
    thispop = pop_prior + dt*r*pop_prior.*(1-pop_prior./K) + h.*pop_prior.*random_numbers.*sqrt(dt) + .5.*pop_prior.*h.^2.*dt.*(random_numbers.^2-1);
    
    if lextinct
        thispop(thispop<popthresh) = 0;
    end
    
    for m = 1:M
     moments(j,m) = mean(thispop.^m);
    end
    
    pop_prior = thispop;
    
    
end

 

end




