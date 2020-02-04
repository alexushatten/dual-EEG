function [Onsets, Offsets] = find_epochs_onsets_offsets( Idx )
% find_epochs_onsets_offsets: find onsets and offsets of epochs defined by
% subscript indices
%
% [Onsets, Offsets] = find_epochs_onsets_offsets( Idx )
%
%  Inputs
% --------
% Idx: subscript indices of epochs
%
%  Outputs
% ---------
% Onsets: subscript indices of epochs onsets
% Offsets: subscript indices of epochs offsets
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, November 2018, validates on test indices 2018-11-27
%-------------------------------------------------------------------------

Diff = dnif(Idx,max(Idx));
Onsets = find(diff(Diff)==1);
Offsets = find(diff(Diff)==-1);
if Diff(end)==1 %numel(Offsets)==(numel(Onsets)-1)
    Offsets(end+1)=length(Diff); % signal cut before the end of epoch
end;
if Diff(1)==1 %numel(Offsets)==(numel(Onsets)+1)
    OnsetsBKP = Onsets;
    Onsets = zeros(length(Onsets)+1,1);
    Onsets(2:end) = OnsetsBKP; % epoch already there when signal starts
end
Onsets = Onsets+1;
if numel(Onsets)~=numel(Offsets)
    error('Onsets & offsets size mismatch')
end

end

