function [fmu,fstd,fCov] = compPredDist_GP(xx,pp)
% [fmu,fstd,fCov] = compPredDist_GP(xx,pp)
%
% INPUT:
% -----
%  xx [nx x d] - spatial input locations in a d-dimensional space (d<=3)
%  pp - parameter struct with fields:
%     .kprs.len - length scale
%     .fdprs.minlen       [1 x d] - minimum length scale for each dimension
%     .fdprs.circinterval [2 x d] - circular support in each stimulus dimension
%     .fdprs.condthresh   [1 x 1] - condition number for thresholding for small eigenvalues 
%     .wmean [m x 1] - mean in Fourier coefficient space
%     .wcov  [m x m] - covariance in Fourier coefficient space
%
% OUTPUT:
% ------
%    fmu [nx x 1]  - mean of predictive distribution points in xx
%   fstd [nx x 1] - 1SD posterior credible intervals for values in ff
%   fCov [nx x nx] - full posterior covariance over ff


% len = pp.kprs.len;
% if isfield(pp.fdprs,'minlen')
%     len = max(len,pp.fdprs.minlen);
% end

% Get Fourier representation for points in xx
FourierBasis = mkRBFfourierBasis(xx,pp.kprs.len,pp.fdprs);

% Compute predictive mean
fmu = FourierBasis*pp.wmean;

% If desired, compute 1SD credible interval
if nargout>1
    fstd = sqrt(sum((FourierBasis*pp.wcov).*FourierBasis,2));
end

% If desired, compute covariance matrix
if nargout>2
    fCov = (FourierBasis*pp.wcov*FourierBasis');
end

