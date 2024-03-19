function fdprs = initFDprs(xvals)
% fdprs = initFDprs(xvals);
%
% Initialize parameters for Fourier-domain representation of GPs
%
% INPUT:
% -------
%    xvals [nt x d] - stimulus matrix (each row contains input point at single time)
%
% OUTPUT:
% -------
%    fdprs (struct) - struct governing Fourier domain representation
%         .minlen [1 x 1] or [m x 1] - minimum length scale for each dimension (can be scalar)
%         .circinterval [2 x d] - circular support in each stimulus dimension
%         .condthresh   [1 x 1] - condition number for thresholding for small eigenvalues 

fdprs.minlen = range(xvals)/10;
fdprs.circinterval = [min(xvals); max(xvals)+5*fdprs.minlen];
fdprs.condthresh = 1e8;

fprintf('-------------------------------------------------------------------------------\n');
fprintf('Initializing Fourier-Domain representation with minimum lengthscale(s) (minlen):\n');
for jj = 1:length(fdprs.minlen)
    fprintf('  %.2f   ', fdprs.minlen(jj));
end
fprintf('\n-------------------------------------------------------------------------------\n');
