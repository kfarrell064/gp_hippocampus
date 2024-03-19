function [ri,ci] = trilinds(nn,k)
%  trilinds(nn,k) - extract row and column indices of lower triangular elements of a matrix
%  of size nn (default k=0 if not provided)
%
%  Inputs:
%   nn - sidelength of square matrix
%    k - which diagonal to start at (0 = main diagonal) (OPTIONAL).
%
%  Outputs:
%   ii - indices of entries of lower triangle (from 1 to nn^2).
%   [ri,ci] - row and column indices of lower triangle

if nargin < 2
    k = 0;
end

if nargout == 1
    ri = find(tril(ones(nn),k));
else
    [ri,ci] = find(tril(ones(nn),k));
end
