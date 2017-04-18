function [params] = mapSigma(params)

params.sigmin = sqrt(2*params.zmin);
params.sigmax = sqrt(2*params.zmax);
params.numsigs = params.numzs;

params.sigmin_an = sqrt(2*params.zmin_an);
params.sigmax_an = sqrt(2*params.zmax_an);
params.numsigs_an = params.numzs_an;



end