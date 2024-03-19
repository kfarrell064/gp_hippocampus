function check_FourierDomainCovariance_nD
% check_FourierDomainCovariance_nD
%
% Check that n-dimensional GP RBF covariance defined in the Fourier domain
% matches that in the real domain to high precision

% Change to parent directory if needed
if contains(pwd,'/unittests')
    cd ..; setpaths; % move up a directory and call setpaths
    addpath unittests;
end
clf;

% 1. Sample random points on the x axis at which to measure the function
nd = 2; % dimensionality of inputs
nx = 200; % number of points
xmax = 10*(1:nd); % maximum of range
xsamp = [sort(rand(nx,1)),rand(nx,nd-1)].*xmax;

% 2. Set GP hyperparameters
rhox = 6.4; % marginal variance
lenx = 2.3; % length scale

% 3. Compute standard GP covariance function
Kx = Krbf(xsamp,lenx,rhox); % x covariance matrix

% 4. Set hyper-parameters governing Fourier-domain representation
xrnge = xmax'*[0 1]; % support of x samples
nSTD = 6; % number of length scales to extend circular interval
fdprs.circinterval = [xrnge(:,1),xrnge(:,2)+nSTD*lenx]'; % circular interval 
fdprs.condthresh = 1e8;  % threshold for cutting off small singular values

%% 5. Compute Fourier-domain covariance
[cxpriordiag,xBasis] = Krbf_fourier(xsamp,lenx,rhox,fdprs); % compute Fourier-domain representation
Kxfd = xBasis*diag(cxpriordiag)*xBasis'; % construct covariance
fprintf('Number of frequencies used: %d\n', length(cxpriordiag));

% -- Make plots showing original and Fourier-domain covariance ----
subplot(221); 
imagesc(Kx); 
title('original covariance'); 
xlabel('x index');ylabel('x index');

subplot(222); 
imagesc(Kxfd);
title('Fourier-domain covariance'); 
xlabel('x index');ylabel('x index');

subplot(223);
iiplot = 1:5:40;
h1 = plot(1:nx,Kx(:,iiplot), 'b'); hold on;
h2 = plot(1:nx, Kxfd(:,iiplot), 'r--'); hold off;
axis tight; box off;
legend([h1(1), h2(1)], 'original', 'Fourier-domain');
title('sample rows'); xlabel('x');

subplot(224);
plot(1:nx,Kx(:,iiplot)-Kxfd(:,iiplot));
title('errors'); xlabel('x');

%% 5. Test if differences are within expected limits

maxdiff = max(abs(Kx(:)-Kxfd(:)));
maxallowedDiff = 10*max(normpdf(nSTD),rhox/fdprs.condthresh);
if maxdiff < maxallowedDiff
        fprintf('\nunit test: PASSED\n');
else
    warning('unit test: FAILED');
end

