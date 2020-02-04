function [Bipoles,ElecMatch,labelsB,labelsB2,EEGb] = bipolar_montage(labels,PairLimit,EEG)
% Make bipolar montage from electrode names and return EEG trace with
% bipolar montage if entered as input (for SEEG, based on channel labels)
% 
% Usage:
%-------------------------------------------------------------------------
% [Bipoles,ElecMatch,labelsB,labelsB2,EEGb] = bipolar_montage(labels,PairLimit,EEG)
%
% [Bipoles,ElecMatch,labelsB,labelsB2] = bipolar_montage(labels,PairLimit)
% 
% Inputs:
%-------------------------------------------------------------------------
%
% labels: cell of strings with electrodes names
%
% PairLimit: maximal "distance" between electrodes such that they form a pair
% (unit is in terms of #electrodes), e.g. should I make a pair with AG1 and
% AG3 given that AG2 is not present in the vector?
%
% EEG: (optional) [channel x time] EEG traces to recalculate based on bipolar
% montage to be created
%
% Outputs:
%-------------------------------------------------------------------------
% 
% Bipoles: index of electrodes to keep
%
% ElecMatch: index of cluster to which each electrode belongs
%
% labelsB: [1 x n] cell of string, labels for bipolar montage, each cell
% contains the original labels concatenated with a "-" (hyphen)
%
% labelsB2: [2 x n] cell of string, labels for bipolar montage but the pair
% of channel is distributed across columns 1 and 2 of the cell array
%
% EEGb: [channel pairs x time] EEG traces recalculated based on bipolar
% montage
%
%-------------------------------------------------------------------------
% NB: not appropriate for scalp EEG bipolar montages!
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, January 2018, updated November 2018
% last update June 2019, fixed issues when some channels are excluded
%-------------------------------------------------------------------------

% if nargin>1
%     EEG = varargin{1};
% end

clear ChNames
for c = 1:length(labels)
    ChNames{c} = char(regexpi(labels{c},'[a-z]+','match')); %#ok<AGROW>
    if size(ChNames{c},1)~=1
        ChNames{c} = ''; %#ok<AGROW>
    end
end
ChNames = ChNames(~cellfun(@isempty,ChNames)); % fix for some cases where channel names are empty or contain several blocks of characters

% here we assumes "Mkr" & empty channels are already out, we just to ignore
% ECG channels:
ChNames = ChNames(cellfun(@isempty,regexpi(ChNames,'^ecg( )*$')))';

ElecMatch = match_vectors(ChNames,unique(ChNames),0);

Bipoles = find(diff(ElecMatch)==0);

clear labelsB
for c = 1:size(Bipoles,1)
    labelsB{c} = [labels{Bipoles(c)},'-',labels{Bipoles(c)+1}]; %#ok<AGROW>
end

clear labelsB2
for c = 1:size(Bipoles,1)
    labelsB2{c,1} = labels{Bipoles(c)}; %#ok<AGROW>
    labelsB2{c,2} = labels{Bipoles(c)+1}; %#ok<AGROW>
end

Keep = true(size(Bipoles,1),1);
if nargin>1
    for c = 1:size(Bipoles,1)
      Keep(c) = ~((str2double(char(regexpi(labelsB2{c,2},'\d','match'))')-str2double(char(regexpi(labelsB2{c,1},'\d','match'))'))>PairLimit);
    end
end
labelsB2 = labelsB2(Keep,:);
labelsB = labelsB(Keep);
Bipoles = Bipoles(Keep);

if nargin>2
%     EEGb = nan(sum(diff(ElecMatch)==0),size(EEG,2));
    EEGb = nan(length(Bipoles),size(EEG,2));
    for c = 1:size(EEGb,1)
        EEGb(c,:) = EEG(Bipoles(c),:)-EEG(Bipoles(c)+1,:);
    end
end

end
