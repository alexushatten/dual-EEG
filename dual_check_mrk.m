function Table = dual_check_mrk(MrkFiles,Fig,ToIgnore)
% dual_check_mrk: get unique instances of markers in Cartool .mrk files,
%   count and make histogram of occurences
%
% Table = dual_check_mrk(MrkFiles,Fig)
%
%  Inputs
% --------
% MrkFiles: [n x 1] cell array with strings of filepaths to "n" .mrk files
% Fig: logical (or numeric), whether or not to plot a figure (default is 0)
% ToIgnore: [n x 1] cell array with strings of markers to ignore (default is '')
%
%  Outputs
% ---------
% Table: markers types with number of occurences, onsets and offsets for
%   each category
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, August 2018
%-------------------------------------------------------------------------

if nargin < 3
    ToIgnore = '';
end

if nargin < 2
    Fig = false;
end

%% Get marker types
EventTypes = {};
for f = 1:numel(MrkFiles)
    [~, ~, MarkerLabel] = read_mrk_file(MrkFiles{f});
    for ff = 1:numel(ToIgnore)
        MarkerLabel = MarkerLabel(cellfun(@isempty,regexp(MarkerLabel,ToIgnore{ff})));
    end
    EventTypes(numel(EventTypes)+1:numel(EventTypes)+numel(unique(MarkerLabel))) = unique(MarkerLabel);
end
EventTypes = unique(EventTypes');

%% Count instances
Counts = zeros(size(EventTypes));
for f = 1:numel(MrkFiles)
    [~, ~, MarkerLabel] = read_mrk_file(MrkFiles{f});
    for ff = 1:numel(ToIgnore)
        MarkerLabel = MarkerLabel(cellfun(@isempty,regexp(MarkerLabel,ToIgnore{ff})));
    end
    Counts = Counts + countmember(EventTypes,MarkerLabel);
end
Table = ([EventTypes,cellstr(num2str(Counts))]);

%% plot figure
if Fig == 1
    %     figure; bar(Counts); set(gca,'xtick',1:numel(Counts)); set(gca,'xticklabel',regexprep(EventTypes,'_',' ')); set(gca,'XTickLabelRotation',90);
    figure; barh(Counts); set(gca,'ytick',1:numel(Counts)); set(gca,'yticklabel',EventTypes,'TickLabelInterpreter','none');
end

%% get indices of occurences:
for f = 1:numel(MrkFiles)
    [MarkerTimeStart, MarkerTimeEnd, MarkerLabel] = read_mrk_file(MrkFiles{f});
    for ff = 1:numel(ToIgnore)
        MarkerLabel = MarkerLabel(cellfun(@isempty,regexp(MarkerLabel,ToIgnore{ff})));
        MarkerTimeStart = MarkerTimeStart(cellfun(@isempty,regexp(MarkerTimeStart,ToIgnore{ff})));
        MarkerTimeEnd = MarkerTimeEnd(cellfun(@isempty,regexp(MarkerTimeEnd,ToIgnore{ff})));
    end
    Idx = match_vectors(EventTypes,MarkerLabel,1);
    for n = 1:size(Idx,1)
        if isa(Idx,'cell')
            Table{n,3}{f} = MarkerTimeStart(Idx{n});
            Table{n,4}{f} = MarkerTimeEnd(Idx{n});
        else
            Table{n,3}{f} = MarkerTimeStart(Idx(n));
            Table{n,4}{f} = MarkerTimeEnd(Idx(n));
        end
    end
end

end