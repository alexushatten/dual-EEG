function z = prc2zscore( p, side )
%% prc2zscore: get z-score necessary to meet percentile of normal distribution (confidence of interval)
%
% z = prc2zscore( p, side )
%
%  Inputs
% --------
% p: [1 x 1] integer, percentile
% side: string, 'unilateral' or 'bilateral' (default: 'bilateral')
%
%  Outputs
% ---------
% z: [1 x 1] integer, z-score corresponding to percentile p of normal
% distribution
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, November 2018
%-------------------------------------------------------------------------

if nargin < 2
    side = 'bilateral';
end
if strcmpi(side,'unilateral')
    p = 1-(1-p)*2;
end
z = erfinv(p)*sqrt(2);

end

