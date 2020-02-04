function RandSampleAepochsNotInB = dual_mrk_rand_equ_xor( Aepochs, Bepochs, PrecisionUnit )
% dual_mrk_rand_equ_xor: extract random samples of continuous data from epochs
% that do not overlap with other epochs with matching length
%
% RandSampleAepochsNotInB = dual_mrk_rand_equ_xor( Aepochs, Bepochs, PrecisionUnit )
%
%  Inputs
% --------
% Aepochs: [n x 1] subscript indices of epochs from which to randomly
%          select
% Bepochs: [n x 1] subscript indices of epochs with which Aepochs should
%          not overlap
% PrecisionUnit: whether to match the length of random selection of epochs
% based on the number of samples ('sample') or epochs ('epoch')
%
%  Outputs
% ---------
% RandSampleAepochsNotInB: random selection of epochs from Aepochs that
% match the length of Bepochs and do not overlap with them
%
%-------------------------------------------------------------------------
% NB: In general, this function works well if there is a lot in Aepochs but
% not so many Bepochs, so we would like to preserve Bepochs as much as
% possible. This would not be for time-locked analysis, though, because
% although epochs of A will be continuous (as much as possible), the output
% of this function results in epochs of variable length and 
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, November 2018
%-------------------------------------------------------------------------

% Exclude periods of Aepochs overlapping with Bepochs:
BnotinA = setdiff(Aepochs,Bepochs);

% Detect onsets and offsets of epochs:
[Onsets, Offsets] = find_epochs_onsets_offsets(BnotinA);

% Random selection:
RandSamples = randsample(1:numel(Onsets),numel(Onsets));
% all(sort(RandSampleRandoms)==unique(RandSampleRandoms)) % there are no repetitions

% Get time frames within Aepochs until we match the number of samples
% for Bepochs:
RandSampleAepochsNotInB = [];
s = 0;
while (length(RandSampleAepochsNotInB) < length(Bepochs)) && (s < length(RandSamples))
    s = s+1;
    RandSampleAepochsNotInB = [RandSampleAepochsNotInB; (Onsets(RandSamples(s)):Offsets(RandSamples(s)))']; %#ok<AGROW>
end
RandSampleAepochsNotInB = sort(RandSampleAepochsNotInB);
if strcmpi(PrecisionUnit,'epoch')
    fprintf('Done. Difference in number of datapoints: %s\n',num2str(length(RandSampleAepochsNotInB)-length(Bepochs)));
elseif strcmpi(PrecisionUnit,'sample')
    if length(RandSampleAepochsNotInB)<length(Bepochs)
        warning('Done, but remaining difference in number of datapoints: A has %s more data points than B !\n',num2str(abs(length(RandSampleAepochsNotInB)-length(Bepochs))));
    else
        RandSampleAepochsNotInB = RandSampleAepochsNotInB(1:length(Bepochs));
    end
else
    error('Unknown precision unit input argument')
end

end

