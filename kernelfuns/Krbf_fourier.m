function [cdiag,fdbasis] = Krbf_fourier(xx,len,rho,fdp)
% Compute basis for RBF covariance in Fourier domain 
%
% [cdiag,fdbasis] = Krbf_fourier(xx,len,rho,fdp)
%
% GP covariance matrix given by: K = Bfft*cdiag*Bfft'
%
% INPUT:
% -----
%           xx [nx x d] - spatial input locations in a d-dimensional space (d<=3)
%           len [1 x 1]  - length scale
%           rho [1 x 1]  - marginal variance
%          fdp (struct) - struct governing Fourier domain representation
%             .circinterval [2 x d] - circular support in each stimulus dimension
%             .condthresh   [1 x 1] - condition number for thresholding for small eigenvalues 
%
% OUTPUT:
% ------
%   cdiag [nb x 1] - vector with thresholded eigenvalues of C
%       B [nb x nx] - column vectors define orthogonal basis for C (on Reals)


% extract size of input data
nd = size(xx,2);  % number of input dimensions

% replicate lens to be same length as # of dimensions
if (nargout==1)
    
    wwsq = mkRBFfourierFreqs(xx,len,fdp); % compute Fourier-domain frequencies only:  (!!! NOTE - this func doesn't exist yet!!!)
    trho = transformRho(rho,len,nd,1); % transformed rho constant
    cdiag =  trho * exp(-.5*(sum(wwsq,2)*len.^2)); % diagonal elements of prior covariance

else

    [fdbasis,wwsq] = mkRBFfourierBasis(xx,len,fdp); % compute Fourier-domain basis and frequencies
    trho = transformRho(rho,len,nd,1); % transformed rho constant
    cdiag =  trho * exp(-.5*(sum(wwsq,2)*len.^2)); % diagonal elements of prior covariance
    
end