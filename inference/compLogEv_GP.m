function logev = compLogEv_GP(xsamp,ysamp,len,rho,nsevar,fdp)
% logev = compLogEv_GP(xsamp,ysamp,len,rho,nsevar,fdp)
%
% Compute log-evidence under GP regression model: 
%
%    log P(y | x, theta) = log \int P(y | f, x) p(f | theta) df

%
% INPUT:
% -------
%       xvals [nt x d] - stimulus matrix (each row contains a single input point)
%       yvals [nt x 1] - response vector (noisy output)
%          len [1 x 1] - length scale 
%          rho [1 x 1] - marginal variance
%       nsevar [1 x 1] - variance of observation noise
%         fdp (struct) - struct governing Fourier domain representation
%            .circinterval [2 x d] - circular support in each stimulus dimension
%            .condthresh   [1 x 1] - condition number for thresholding for small eigenvalues 
%
% OUTPUT:
% --------
%   logev - log evidence or marginal likelihood: 

nsamp = size(xsamp,1);
% Compute posterior using Fourier-domain version
[cdiag,fBasis] = Krbf_fourier(xsamp,len,rho,fdp);  % diagonal prior and basis

% Fourier-domain formula
XX = fBasis'*fBasis/nsevar;
XY = fBasis'*ysamp/nsevar;
Cinv = diag(1./cdiag);
trm1 = -.5*(logdet(XX+Cinv) + sum(log(cdiag)) + nsamp*log(2*pi*nsevar)); % Log-determinant term
trm2 = .5*(-ysamp'*ysamp/nsevar + XY'*((XX+Cinv)\XY));   % Quadratic term
logev = trm1+trm2;  % negative log evidence
