% lesExactStationaryMoments.m
%
% Implement analytic results for LES stationary moments w/out extinction,
% see Schlomann (2017).
%
% Date:  4/10/17 - First written


function [momentsAn] = lesExactStationaryMoments(r,K,sig,M)

% pre-allocate
momentsAn = zeros(1,M);

% Compute each moment individually, for clarity
for m = 1:M
    
    % for product term
    mvec = 1:(m-1);
    prodfactor = 1+(mvec-1).*sig^2./r./2;
    
    % the mth moment
    % For large m, K, can remove K^m to compute moments of X/K, avoiding large
    % numbers
    momentsAn(m) = K^m.*(1-sig.^2/r./2).*prod(prodfactor);

end