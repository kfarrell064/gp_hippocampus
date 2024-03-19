function K = Krbf_circular(x,len,rho,xperiod)
% Make covariance matrix using 1-dimensional RBF covariance function (aka
% squared-exponential kernel) on a circular boundary by summing up 
% distances between points in circular space
%
% K = Krbf_cirrcular(x,prs)
%
% Covariance matrix parametrized as:
%  C_ij = rho*exp(((x_i-x_j)^2/(2*len^2))  
%
% INPUTS:
%     x [nx x 1] - points at which to evaluate kernel function (assumes all
%                  live within a single period)
%     len - length scale
%     rho - marginal variance
%   xcirc - circular interval
%
% OUTPUT:
%   K [nx x nx] - covariance matrix

% Initialize covariance
K = zeros(length(x));

% determine # periods to use based on keeping dists <= 3 lengthscales
nperiods = 1+ceil(3*len/xperiod);

if nperiods > 50
    warning('Krbf_circular: lengthscale >> period');
    fprintf('\nrank(K) approaching 1!\n');
    fprintf('\n(Using nperiods=50 to compute K)\n');
    nperiods  = 50;
end

% loop over periods
for jj = -nperiods:nperiods
    % Compute squared distances for this particular shift of the data
    sqrdists = sum(bsxfun(@minus,(x(:)+jj*xperiod),x(:)').^2,3);

    % Compute contribution to covariance 
    K = K + rho*exp(-.5*sqrdists/len.^2); % the covariance matrix
end
    
