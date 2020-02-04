function Idx = match2bins_before( V, Bins )
% match2bins: get indices of where each element of V is mapped to Bins,
% based on smallest absolute difference, useful for mapping values to bins
% => this version of the function forces to match the bin that is smaller
% than each value
%
% Idx = match2bins( V, Bins )
%
%  Inputs
% --------
% V: [n x 1] vector of values to map
% Bins: [n x 1] vector of values to map to
%
%  Outputs
% ---------
% Idx: indices of V in Bins
%
%-------------------------------------------------------------------------
% NB: equivalent of match_vectors.m for non-exact solutions
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, May 2019
%-------------------------------------------------------------------------

Idx = nan(size(V));
for i = 1:length(V)
    Temp = V(i)-Bins;
    Temp(Temp<0)=[];
    [~,IdxTemp] = min(Temp);
    Idx(i) = IdxTemp;
end

end

