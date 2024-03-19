function check_FourierDomainCovariance_ungridded
% check_FourierDomainCovariance_ungridded
%
% Check that 1D GP RBF covariance defined in the Fourier domain matches 
% that in the real domain to high precision, for ungridded x values

% Change to parent directory if needed
if contains(pwd,'/unittests')
    cd ..; setpaths; % move up a directory and call setpaths
    addpath unittests;
end
clf;

% 1. Set GP hyperparameters
rhox = 11.3; % marginal variance
lenx = 3.2; % length scale

% 2. Sample random points on the x axis at which to measure the function
nx = 200; % number of points
xmax = 20; % maximum of range
xsamp = sort(rand(nx,1)*xmax); 

% 3. Compute standard GP covariance function
Kx = Krbf(xsamp,lenx,rhox); % x covariance matrix

% 4. Set hyper-parameters governing Fourier-domain representation
xrnge = [0,xmax]; % support of x samples
nSTD = 6; % number of length scales to extend circular interval
fdprs.circinterval = [xrnge(1),xrnge(2)+nSTD*lenx]'; % circular interval 
fdprs.condthresh = 1e8;  % threshold for cutting off small singular values

% 5. Compute Fourier-domain covariance
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
h1 = plot(xsamp,Kx(:,iiplot), 'b'); hold on;
h2 = plot(xsamp, Kxfd(:,iiplot), 'r--'); hold off;
axis tight; box off;
legend([h1(1), h2(1)], 'original', 'Fourier-domain');
title('sample rows'); xlabel('x');

subplot(224);
plot(xsamp,Kx(:,iiplot)-Kxfd(:,iiplot));
title('errors'); xlabel('x');

%% 5. Test if differences are within expected limits

maxdiff = max(abs(Kx(:)-Kxfd(:)));
maxallowedDiff = 10*max(normpdf(nSTD),rhox/fdprs.condthresh);
if maxdiff < maxallowedDiff
        fprintf('\nunit test: PASSED\n');
else
    warning('unit test: FAILED');
end

