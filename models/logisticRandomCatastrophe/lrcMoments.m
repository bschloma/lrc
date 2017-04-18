
% Program:  lrcMoments.m
% 
% Summary:  Compute moments of Logistic Random Catastrophe model.  Populations grow logistically but are
%           subject to Poisson catastrophes of constant fraction.  No
%           extinction.
% 
% Inputs
%           mu      : Bacterial growth rate (default 1.5 hr^{-1})
%           Kparams : Bacterial carrying capacity .
%                     Can be a single number, or [mean K, sigma of *log10*K for a
%                     log-normal distribution]. In the latter case, each run will
%                     draw K from a log-normal distribution. Default [1e4 0.5]
%           f       : post-collapse population fraction
%           lambda  : probability per unit time of population collapse (default
%                     0.1 hr^{-1})
%           numtrials  : number of simulation runs
%           Tmax    : hrs., end time of simulation (default 24 hours)
%           dt      : hrs., time step of simulation (default 0.1 hours)
%           lextinct: Logical for extinction, true=popthresh (fixed at 1)
%                     is an absorbing state, false = no extinction
%           M       : Integer, number of moments (>0)
% 
% 
% Outputs
%           moments   : Ensemble moments; row = time; columns = moment #; 
%           
% 
% Fixed parameters
%           Initial population size : 10
% 
% Author:   Brandon Schlomann, modified from stochastic_collapse_model by Raghuveer Parthasarathy
% 
% Date:     3/18/17

function [moments] = lrcMoments(mu,Kparams,f, lambda, numtrials, Tmax, dt,lextinct,M)

rng('shuffle'); %initialize random number generator

% Default parameter values
if ~exist('lambda', 'var') || isempty(lambda)
    lambda = 0.1; % hr^{-1}, probability per unit time of collapse
end
if ~exist('mu', 'var') || isempty(mu)
    mu = 1.5; % 1/hours
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
    numtrials = 1000;
end

if ~exist('f', 'var') || isempty(f)
    f = .01;
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
pop0 = sqrt(Kparams(1));%10;

pop_prior = pop0*ones(1,numtrials);
moments(1,:) = pop0.^mvec;

% The actual simulation
for j=2:numsteps
     
    random_numbers = rand(1,numtrials); % pre-calculate random numbers, [0, 1]
    dNt = random_numbers < (lambda*dt); % true at time points where a collapse occurs
    
    thispop = (1-dNt).*pop_prior + (1-dNt).*dt*mu.*pop_prior.*(1-pop_prior./K) + dNt.*f.*pop_prior;
    
    if lextinct
        thispop(thispop<popthresh) = 0;
    end
    
    for m = 1:M
     moments(j,m) = mean(thispop.^m);
    end
    
    pop_prior = thispop;
    
    
end

 

end




