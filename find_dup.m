function Dup = find_dup( V )
% find_dup: find duplicates after the first occurence in a vector
%
% Dup = find_dup( V )
%
%  Inputs
% --------
% V : 1-dimensional vector
%
%  Outputs
% ---------
% Dup : indices of duplicates after first occurence
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, June 2018
%-------------------------------------------------------------------------

Ndim = ndims(V);
if Ndim ~= 2 || ~any(size(V)==1)
    error('Function not defined for arrays or if number of dimensions is not 2')
end

V = V(:);

[~, I] = unique(V, 'first');
Dup = 1:length(V);
Dup(I) = [];

end

