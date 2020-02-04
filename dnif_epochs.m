function xdi = dnif_epochs(idxStart, idxEnd,L)
% DNIF_EPOCHS (aka reverse find): get logical vector from subscript indices
% obtained with find.m, similar to dnif.m but for epochs
%
%   idxStart:    vector of subscript indices for epochs onset
%   idxEnd:      vector of subscript indices for epochs offset
%   L:      length of the logical vector requested
%
% Renaud Marquis @ FBM lab, November 2018

if length(idxStart)~=length(idxEnd)
    error('Onsets / offsets length mismatch')
end

xdi = false(L,1);

for b = 1:length(idxStart)
    xdi(idxStart(b):idxEnd(b)) = true;
end

end