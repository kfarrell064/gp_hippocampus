function [neglogev,grad,H,mupost,Lpost,ii,Cinv] = neglogEv_GPregress(prs,dd,wwsq,condthresh)
% [neglogev,grad,H,mupost,Lpost,ii,Cinv] = neglogEv_GPregress(prs,dd,wwsq,condthresh)
%
% Computes negative log-evidence: 
%
%    -log P(Y|length, rho, sig^2) 
%
% under linear-Gaussian model: 
%
%       f ~ N(0,C),
%       y = f + n,  n ~ N(0,sig^2)
%
% where C is RBF or "squared exponential" covariance matrix
% 
% INPUT:
% -------
%        prs [3 x 1] - RBF-GP parameters [log(len); log(rho); log(nsevar)]
%        dd [struct] - sufficient statistics for regression:
%                      .xx - Fourier domain stimulus autocovariance matrix X'*X 
%                      .xy - Fourier domain stimulus-response cross-covariance X'*Y 
%                      .yy - response variance Y'*Y
%                      .ny - number of samples 
%       wwsq [m x 1] - vector of normalized DFT frequencies, along each dimension
% condthresh [1 x 1] - threshold for condition number of K (Default = 1e8).
%
% OUTPUT:
% -------
%   neglogev - negative marginal likelihood
%   grad - gradient
%   H - Hessian
%   mupost - mean of posterior over regression weights
%   Lpost - posterior covariance over regression weights
%   ii - logical vector indicating which DFT frequencies are not pruned
%   Cinv - inverse prior covariance in diagonalized, pruned Fourier space


% Unpack parameters
len = exp(prs(1));
trho = exp(prs(2));
nsevar = exp(prs(3));

% Compute diagonal representation of prior covariance (Fourier domain)
ii = wwsq < 2*log(condthresh)/len^2;
ni = sum(ii); % rank of covariance after pruning

% Build prior covariance matrix from parameters
% (Compute derivatives if gradient and Hessian are requested)
switch nargout 
    case {0,1} % compute just diagonal of C
        cdiag = Krbf_fourierdiag(len,trho,wwsq(ii)); % compute diagonal and Fourier freqs

    case 2
        % compute diagonal of C and deriv of C^-1
        [cdiag,dcinv] = Krbf_fourierdiag(len,trho,wwsq(ii)); % compute diagonal and Fourier freqs

    otherwise
        % compute diagonal of C, 1st and 2nd derivs of C^-1 and C
        [cdiag,dcinv,dc,ddcinv] = Krbf_fourierdiag(len,trho,wwsq(ii)); % compute diagonal and Fourier freqs
end


% Prune XX and XY Fourier coefficients and divide by nsevar
Cinv = spdiags(1./cdiag,0,ni,ni); % inverse cov in diagonalized space
XX = dd.xx(ii,ii)/nsevar; 
XY = dd.xy(ii)/nsevar;

% Compute negative log-evidence
trm1 = -.5*(logdet(XX+Cinv) + sum(log(cdiag)) + (dd.nsamps)*log(2*pi*nsevar)); % Log-determinant term
trm2 = .5*(-dd.yy/nsevar + XY'*((XX+Cinv)\XY));   % Quadratic term
neglogev = -trm1-trm2;  % negative log evidence

% Compute neglogli and Gradient
if nargout >= 2  
    % Compute matrices we need
    Lpostinv = (XX+Cinv);
    Lpost = inv(Lpostinv);
    Lpdiag = diag(Lpost);
    mupost = Lpostinv\XY;
        
    % --- Compute gradient ------------
    % Derivs w.r.t hyperparams rho and len
    dLdthet = -.5*dcinv'*(cdiag - (Lpdiag + mupost.^2));
    % Deriv w.r.t noise variance 'nsevar'
    RR = .5*(dd.yy/nsevar - 2*mupost'*XY + mupost'*XX*mupost)/nsevar; % Squared Residuals / 2*nsevar^2
    Tracetrm = .5*(ni-dd.nsamps-sum((Lpdiag./cdiag)))/nsevar;
    dLdnsevar = -Tracetrm-RR;
    % Combine them into gardient vector
    expprs = [len;trho;nsevar];
    grad = [dLdthet; dLdnsevar].*expprs;
end

% Compute Hessian 
if nargout >= 3 

    % theta terms (rho and len)
    nthet = 2;  % number of theta variables (rho and len)
    vn = ones(1,nthet); % vector of 1s of length ntheta
    dLpdiag = -(Lpost.^2)*dcinv; % Deriv of diag(Lpost) w.r.t thetas
    dmupost = -(Lpost)*bsxfun(@times,dcinv,mupost); % Deriv of mupost w.r.t thetas
    [ri,ci] = triuinds(nthet);  % get indices for rows and columns of upper triangle
    trm1stuff = -.5*(dc - (dLpdiag + 2*dmupost.*mupost(:,vn)));
    ddLddthet_trm1 = sum(trm1stuff(:,ri).*dcinv(:,ci),1)';
    ddLddthet_trm2 = -.5*ddcinv'*(cdiag - (Lpdiag + mupost.^2));
    ddLddthet = ddLddthet_trm1+ddLddthet_trm2;

    % nsevar term
    dLpdiagv = sum(Lpost.*(Lpost*XX),2)/nsevar; % Deriv of diag(Lpost) wr.t nsevar
    dmupostv = -(Lpost*(mupost./cdiag))/nsevar; % Deriv of mupost w.r.t nsevar
    ddLdv = -(dLdnsevar/nsevar - RR/nsevar ...
        - sum(dLpdiagv./cdiag)/(2*nsevar) ...
        +((-XY+XX*mupost)'*dmupostv)/nsevar);  % 2nd deriv w.r.t. nsevar
        
    % Cross term theta - nsevar
    ddLdthetav = .5*dcinv'*(dLpdiagv+2*dmupostv.*mupost);
    
    % assemble Hessian 
    Hprs = unvecSymMtxFromTriu([ddLddthet;ddLdthetav; ddLdv]);
    H = (Hprs.*(expprs*expprs') + diag(grad));
    
end
