function check_FourierDomainCovariance_gridded
% function check_FourierDomainCovariance_gridded
%
% Check that 1D GP RBF covariance defined in the Fourier domain matches
% that in the real domain to high precision, for gridded x values

% Change to parent directory if needed
if contains(pwd,'/unittests')
    cd ..; setpaths; % move up a directory and call setpaths
    addpath unittests;
end
clf;

% 1. Set GP hyperparameters
rhox = 13.78; % marginal variance
lenx = 8.4; % length scale

% 2. Set a regular grid of x points at which to measure the function
nx = 200; % number of points
xsamp = (1:nx)'; % regular grid 

% 3. Compute standard GP covariance function
Kx = Krbf(xsamp,lenx,rhox); % x covariance matrix

% 4. Set hyperparameters governing Fourier-domain covariance
xrnge = [0,nx]'; % support of x samples
nSTD = 6; % number of length scales to extend circular interval
fdprs.circinterval = [xrnge(1);xrnge(2)+nSTD*lenx]; % circular interval 
fdprs.condthresh = 1e8;  % threshold for cutting off small singular values

% 5. Compute Fourier-domain covariance
[cxpriordiag,xBasis] = Krbf_fourier(xsamp,lenx,rhox,fdprs); % compute Fourier-domain version
Kxfd = xBasis*diag(cxpriordiag)*xBasis';
fprintf('Number of frequencies used: %d\n', length(cxpriordiag));

% -- Make plots showing original and Fourier-domain covariance ----
subplot(221); 
imagesc(Kx); xlabel('x');ylabel('x');
title('original covariance'); axis square;

subplot(222); 
imagesc(Kxfd); xlabel('x');ylabel('x');
title('Fourier-domain covariance'); axis square;

subplot(223);
iiplot = 1:5:40; % columns to plot
h1 = plot(1:nx,Kx(:,iiplot), 'b'); hold on;
h2 = plot(1:nx, Kxfd(:,iiplot), 'r--'); hold off;
axis tight; box off;
legend([h1(1), h2(1)], 'original', 'Fourier-domain');
title('sample rows'); xlabel('x');

subplot(224);
plot(1:nx,Kx(:,iiplot)-Kxfd(:,iiplot));
title('errors'); xlabel('x');

%% 6. Test if differences are within expected limits

maxdiff = max(abs(Kx(:)-Kxfd(:)));
maxallowedDiff = 10*max(normpdf(nSTD),rhox/fdprs.condthresh);
if maxdiff < maxallowedDiff
        fprintf('\nunit test: PASSED\n');
else
    warning('unit test: FAILED');
end

