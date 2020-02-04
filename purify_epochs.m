function [ MOIstart, MOIend, MOIlabel ] = purify_epochs( MarkerStart, MarkerEnd, MarkerLabel, PatternsToKeep, PatternsToExclude, TimeDiffThresh, Direction )
% purify_epochs: from markers of interest (MOI) and markers to exclude,
% filter MOIs to keep only those with no overlap with MTEs using time
% difference threshold of TimeDiffThresh
%
% [ MOIstart, MOIend, MOIlabel ] = ...
%       purify_epochs( MarkerStart, MarkerEnd, MarkerLabel, ...
%           PatternsToKeep, PatternsToExclude, TimeDiffThresh )
%
%  Inputs
% --------
% MarkerStart      : [n x 1] timings of markers' start
% MarkerEnd        : [n x 1] timings of markers' end
% MarkerLabel      : {n x 1} markers labels
% PatternsToKeep    : cell of strings corresponding to markers labels to
%                     keep, regular expressions are supported
% PatternsToExclude : cell of strings corresponding to markers labels to 
%                     exclude, regular expressions are supported
% TimeDiffThresh    : integer, time difference (same units as MarkersStart
%                     & MarkersEnd !) threshold under which markers are
%                     considered as overlapping, i.e. any triplet of MarkerStart,
%                     MarkerEnd & MarkerLabel matching PatternsToKeep will
%                     be kept unless the time difference between them and
%                     another triplet of MarkerStart, MarkerEnd & MarkerLabel
%                     that matches PatternsToExclude is below TimeDiffThresh.
% Direction         : string, whether to exclude markers for which there is 
%                     an overlapping marker after ('forward'), before
%                     ('backward') or both before and after ('both')
%
%  Outputs
% ---------
% MOIstart       : [n x 1] timings of start of filtered markers
% MOIend         : [n x 1] timings of end of filtered markers
% MOIlabel       : {n x 1} labels of filtered markers
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, November 2018, validated on test marker samples 2018-11-27
%-------------------------------------------------------------------------

MOIstart = MarkerStart;
MOIend = MarkerEnd;
MOIlabel = MarkerLabel;

% marker indices space
KeepThem = [];
for p = 1:numel(PatternsToKeep)
    KeepThem = [KeepThem; find(~cellfun(@isempty,regexp(MOIlabel,PatternsToKeep{p})))]; %#ok<AGROW>
end
KeepThem = sort(KeepThem);
if isempty(KeepThem)
    error('No marker matching any pattern was found !')
end

ExcludeThem = [];
for p = 1:numel(PatternsToExclude)
    ExcludeThem = [ExcludeThem; find(~cellfun(@isempty,regexp(MOIlabel,PatternsToExclude{p})))]; %#ok<AGROW>
end
ExcludeThem = sort(ExcludeThem);

% timings space
TimeStartToExclude = MOIstart(ExcludeThem);
TimeEndToExclude = MOIend(ExcludeThem);
MaxTime = max(MOIend)+TimeDiffThresh;
if strcmpi(Direction,'both')
    TimeEndToExclude = TimeEndToExclude+TimeDiffThresh;
    TimeStartToExclude = TimeStartToExclude-TimeDiffThresh;
elseif strcmpi(Direction,'forward')
    TimeEndToExclude = TimeEndToExclude+TimeDiffThresh;
elseif strcmpi(Direction,'backward')
    TimeStartToExclude = TimeStartToExclude-TimeDiffThresh;
end
TimeEndToExclude(TimeStartToExclude<1)=[];
TimeStartToExclude(TimeStartToExclude<1)=[];

TimeStartToKeep = MOIstart(KeepThem);
TimeEndToKeep = MOIend(KeepThem);
LabelsToKeep = MOIlabel(KeepThem);
OnOff2keep = dnif_epochs(TimeStartToKeep,TimeEndToKeep,MaxTime);
OnOff2exclude = dnif_epochs(TimeStartToExclude,TimeEndToExclude,MaxTime);
Overlap = (OnOff2keep + OnOff2exclude) > 1;
OnOff2keep(Overlap) = false;
ToKeep = find(OnOff2keep);
ToDiscard = setdiff(TimeStartToKeep,ToKeep);

DiscardThem = nan(length(ToDiscard),1); % normally should not exceed this size, but sometimes there are two markers with the same start time (and looking and end time does not necessarily help), so we will increase the array size in case this happens...
Count = 0;
for t = 1:length(ToDiscard)
    Count = Count+1;
    DiscardThem(Count:Count+length(find(TimeStartToKeep==ToDiscard(t)))-1) = find(TimeStartToKeep==ToDiscard(t));        
end

MOIstart = TimeStartToKeep;
MOIend = TimeEndToKeep;
MOIlabel = LabelsToKeep;

MOIstart(DiscardThem)=[];
MOIend(DiscardThem)=[];
MOIlabel(DiscardThem)=[];

if isempty(MOIstart) || isempty(MOIend) || isempty(MOIlabel)
    warning('Markers matching patterns were found, but after exclusion of markers overlapping with markers matching patterns to exclude, none of them remain...')
end

end

