function mkr = dual_clean_mrk_NEW(markers1,markers2,MustOverlap,Period,DINtype)

% for NEW REALIGN, 2018-12-05
% GOAL: 
%   - quality check on marker files for icEEG (markers1) and hdEEG (markers 2)
%
% NOTE:
%   - use dual_read_mkr to read textfile into matlab
%   - labels from micromed only go up to 256 and then start again at 1, LS
%   has relabeled markers above 256.
%
% INPUT: marker filenames for icEEG and HD-EEG
%        MustOverlap: in some cases, the marker labels are not consistent
%        but the files do overlap => in this case, the issue is that marker
%        labels should be unwrapped for one of the marker files because it
%        is already unwrapped for the other...
%
%        If specifying MustOverlap, Period must be defined (e.g. 255)
%
%        DINtype: 'DINin' or 'DINafter' (whether the number to extract is "128" or
%        "12" in e.g. "DIN12:128")
%
% OUTPUT: common marker file 
%   - column 1: marker label
%   - column 2: icEEG datapoints
%   - column 3: hdEEG datapoints
%
% DEVELOPMENT
%   - MArch 2017 Maxime Baud
%   - Aug 2017, Renaud Marquis, last updated July 2019

%#RM, 04.07.17 18:10:48
folder1 = fileparts(markers1);
folder2 = fileparts(markers2);
if ~strcmp(folder1,folder2)
    warning('hdEEG and icEEG not in the same folder, log file will be saved in folder containing icEEG file, i.e. %s',folder1)
end
if isempty(folder1)
    folder = cd;
else
    folder = folder1;
end
fid = 0; % dirty fix, 2018-12-05
% fid = fopen(fullfile(folder,'Marker_file_check_summary_log.txt'),'w'); % print log of what was realigned
% fprintf(fid,'Created: %s \n\n',datestr(now,'dd-mmm-yyyy'));
% 
% if ~strcmp(folder1,folder2)
%     fprintf(fid,'hdEEG and icEEG not in the same folder, log file will be saved at icEEG file folder, i.e. %s\n\n',folder1);
%     warning('hdEEG and icEEG not in the same folder, log file will be saved at icEEG file folder, i.e. %s',folder1);
% end
% fprintf(fid,'Full paths: \nicEEG file: %s \nhdEEG file: %s\n\n',markers1,markers2);

WarningMessage = zeros(2,10);

%% PART 1 ====== CLEAN LABELS

% check for non-numeric markers, multiple ones, doublets, gaps, extra
% markers (almost but not doublets), markers only in one file, and frames
% with sampling rate issues.

[mkr, NonNumMark, multi_ones, MultiOnes, delta, DoubletsGaps, ExtraMarkersDiscarded, icEEGbutNOThdEEG, hdEEGbutNOTicEEG, mkr1, mkr2, idx, markers1,markers2,folder1,folder2,folder,fid,WarningMessage] = clean_mkr_part1(markers1,markers2,folder1,folder2,folder,fid,WarningMessage,DINtype); %#ok<ASGLU>

if MustOverlap
    warning('Assuming files should overlap based on period of %d, will wrap marker labels if one of the files is unwrapped and the other not !')
    [mkr, NonNumMark, multi_ones, MultiOnes, delta, DoubletsGaps, ExtraMarkersDiscarded, icEEGbutNOThdEEG, hdEEGbutNOTicEEG, mkr1, mkr2, idx, markers1,markers2,folder1,folder2,folder,fid,WarningMessage] = clean_mkr_part1_but_deunwrap_before(markers1,markers2,folder1,folder2,folder,fid,WarningMessage,Period,DINtype); %#ok<ASGLU>
end

%% PART 2 ====== LABELS UNWRAPPING
% At this stage, normally, there should not be any duplicate markers
% anymore!
% Time to check if marker labels were unwrapped or not beforehand...
if any(diff(mkr1(:,2))<0) || any(diff(mkr2(:,2))<0)
    warning('...UNWRAPPING MARKER LABELS...')
    
    if ~MustOverlap
        % START AGAIN FROM 0 (because removed missing markers should have been kept)
        % and UNWRAP MARKER LABELS:
        [mkr, delta, mkr1, mkr2, DoubletsGaps, ExtraMarkersDiscarded, markers1,markers2,folder1,folder2,folder,fid,WarningMessage] = clean_mkr_unwrap(markers1,markers2,folder1,folder2,folder,fid,WarningMessage,DINtype); %#ok<ASGLU>
    else
        [mkr, delta, mkr1, mkr2, DoubletsGaps, ExtraMarkersDiscarded, markers1,markers2,folder1,folder2,folder,fid,WarningMessage] = clean_mkr_deunwrap_and_unwrap(markers1,markers2,folder1,folder2,folder,fid,WarningMessage,Period,DINtype); %#ok<ASGLU>
        
    end
    
    % THEN CLEAN LABELS AGAIN:
    [mkr, multi_ones, MultiOnes, icEEGbutNOThdEEG, hdEEGbutNOTicEEG, mkr1, mkr2, idx, markers1,markers2,folder1,folder2,folder,fid,WarningMessage] = clean_mkr_part2(markers1,markers2,folder1,folder2,folder,fid,WarningMessage,mkr1,mkr2); %#ok<ASGLU>
    
end

% Merge markers
mkr = [mkr1(:,2), mkr1(:,1), mkr2(:,1)];

%% PART 3 ====== CLEAN DATAPOINTS

% if ~(size(mkr,1)<2)
% C = diff(mkr);
% fs1 = median(C(:,2));
% fs2 = median(C(:,3));
% 
% % Period mismatch
% % period between icEEG and hdEEG is not the same
% % dc = abs(C(:,2)/fs1)./C(:,1)-(C(:,3)/fs2)>0.0020; %#RM, Aug 2017: in this
% % case, means that only 1 sample more or less is allowed on both sides maximum!
% SRthreshold = 0.001; % only 1 millisecond more (or less) on maximum both sides
% %#RM, 2017-08-25: other method:
% dc1 = ((abs(C(:,2)/fs1)./C(:,1)-abs(C(:,3)/fs2)./C(:,1))>SRthreshold*2); %#RM, 2017-08-25
% dc2 = ((abs(C(:,2)/fs1)./C(:,1)-abs(C(:,3)/fs2)./C(:,1))<-SRthreshold*2); %#RM, 2017-08-25
% dc = (dc1+dc2)>0;
% % #RM-2017-08-25: add an exception: sometimes one marker is delayed and then by summing the
% % two frames around this marker, the total duration equals sampling rate x
% % 2 => in this case, shift the marker in time such that it cuts the "double
% % period with shifted marker" in two (~)equals parts (to one sample)
% DC = dc1-dc2;
% DCok = find(((diff(DC)>1)+(diff(DC)<(-1)))>0);
% 
% for n = 1:length(DCok)
%     TwoPeriodsTest = sum([C(DCok(n):DCok(n)+1,2)/fs1 C(DCok(n):DCok(n)+1,3)/fs2])./sum(C(DCok(n):DCok(n)+1,1));
%     if all(abs(TwoPeriodsTest-1)<SRthreshold) % same constraint as above (1 sample on any direction on each side max)
%         [~,IdxTPT] = max(abs(TwoPeriodsTest-1));
%         mkr(DCok(n)+1,IdxTPT+1)=round(mkr(DCok(n),IdxTPT+1)+(mkr(DCok(n)+2,IdxTPT+1)-mkr(DCok(n),IdxTPT+1))/sum(C(DCok(n):DCok(n)+1,1)));
%         dc(DCok(n):DCok(n)+1)=0;
%     end
% end
% 
% if sum(dc)~=0
%     WarningMessage(1:2,8)=1;
%     warning('Some frames had different length in the 2 files and were discarded')
%     disp([mkr(dc,1) C(dc,2:3)])
%     mkr(dc,:)=[];
% end
% PeriodMismatch = find(dc);
% 
% % Long period - allow 3xfs additional periods
% % Probably occasional delay in trigger send when passes maximum of 255
% % pass = (C(:,2)>fs1-4 & C(:,2)<fs1+4) & (C(:,3)>fs2-3 & C(:,3)<fs2+3);
% C = diff(mkr); %#RM Aug 2017: necessary to re-do it!
% dc = abs(C(:,2)-fs1) > 3*fs1 | abs(C(:,3)-fs2) > 3*fs2; %#RM: 3 or 4 seconds?
% if sum(dc)~=0
%     WarningMessage(1:2,9)=1;
%     warning('Some frames are longer than the sampling frequency, discarded frames > 4 sec')
%     displ= abs(C(:,2)-fs1) > fs1 | abs(C(:,3)-fs2) > fs2;
%     disp([mkr(displ,1) C(displ,2:3)])
%     mkr(dc,:)=[];
% end
% LongPeriods = find(dc);
% else
%     PeriodMismatch = [];
%     LongPeriods = [];
% end
% % if logical(FlagRemovedOneSetFromICEEG)
% %     WarningMessage(1,10)=1;
% % elseif logical(FlagRemovedOneSetFromHDEEG)
% %     WarningMessage(2,10)=1;
% % end
% 
% WarningMessage = WarningMessage';
% Strings = {'Discarding non-numerical markers'
%     'Discarding multiple 1 markers'
%     'Discarding marker doublets'
%     'Gaps in the markers'
%     'Discarding extra markers'
%     'Discarding markers from file 1 not in file 2'
%     'Discarding markers from file 2 not in file 1'
%     'Frames with unequal length in 2 files'
%     'Frames longer than sampling rate'
%     'Label unwrapping'
%     };

% fprintf(fid,'\tError types \t\t\t\t\t\t\t\t\ticEEG file \thdEEG file \t#cases icEEG \t#cases hdEEG\n');
% fprintf(fid,['\t',Strings{1},'\t\t\t\t',num2str(WarningMessage(1,1)),'\t\t\t',num2str(WarningMessage(1,2)),'\t\t\t',num2str(length(NonNumMark{1})),'\t\t\t\t',num2str(length(NonNumMark{2})),'\n']);
% fprintf(fid,['\t',Strings{2},'\t\t\t\t\t',num2str(WarningMessage(2,1)),'\t\t\t',num2str(WarningMessage(2,2)),'\t\t\t',num2str(length(MultiOnes{1})),'\t\t\t\t',num2str(length(MultiOnes{2})),'\n']);
% fprintf(fid,['\t',Strings{3},'\t\t\t\t\t\t',num2str(WarningMessage(3,1)),'\t\t\t',num2str(WarningMessage(3,2)),'\t\t\t',num2str(sum(DoubletsGaps{1}==0)),'\t\t\t\t',num2str(sum(DoubletsGaps{2}==0)),'\n']);
% fprintf(fid,['\t',Strings{4},'\t\t\t\t\t\t\t\t',num2str(WarningMessage(4,1)),'\t\t\t',num2str(WarningMessage(4,2)),'\t\t\t',num2str(sum(ExtraMarkersDiscarded{1}>1)),'\t\t\t\t',num2str(sum(ExtraMarkersDiscarded{2}>1)),'\n']);
% % WarningMessage(5,:) removed, redundant with DoubletGaps, #RM-2017-08-21
% fprintf(fid,['\t',Strings{6},'\t',num2str(WarningMessage(6,1)),'\t\t\t',num2str(WarningMessage(6,2)),'\t\t\t',num2str(length(icEEGbutNOThdEEG)),'\t\t\t\t','\n']);
% fprintf(fid,['\t',Strings{7},'\t',num2str(WarningMessage(7,1)),'\t\t\t',num2str(WarningMessage(7,2)),'\t\t\t','\t\t\t\t',num2str(length(hdEEGbutNOTicEEG)),'\n']);
% fprintf(fid,['\t',Strings{8},'\t\t\t',num2str(WarningMessage(8,1)),'\t\t\t',num2str(WarningMessage(8,2)),'\t\t\t',num2str(length(PeriodMismatch)),'\t\t\t\t',num2str(length(PeriodMismatch)),'\n']);
% fprintf(fid,['\t',Strings{9},'\t\t\t\t',num2str(WarningMessage(9,1)),'\t\t\t',num2str(WarningMessage(9,2)),'\t\t\t',num2str(length(LongPeriods)),'\t\t\t\t',num2str(length(LongPeriods)),'\n']);
% fprintf(fid,['\t',Strings{10},'\t\t\t\t\t\t\t\t',num2str(WarningMessage(10,1)),'\t\t\t',num2str(WarningMessage(10,2)),'\t\t\t',num2str(sum(ExtraMarkersDiscarded{1}<0)),'\t\t\t\t',num2str(sum(ExtraMarkersDiscarded{2}<0)),'\n']);
% 
% fclose(fid);

% save(fullfile(folder,'Marker_file_check_summary_log.mat'),...
%     'NonNumMark','MultiOnes','DoubletsGaps','ExtraMarkersDiscarded',...
%     'icEEGbutNOThdEEG','hdEEGbutNOTicEEG','PeriodMismatch','LongPeriods')

end

%% /////// INTERNAL FUNCTIONS ///////

function [mkr, NonNumMark, multi_ones, MultiOnes, delta, DoubletsGaps, ExtraMarkersDiscarded, icEEGbutNOThdEEG, hdEEGbutNOTicEEG, mkr1, mkr2, idx, markers1,markers2,folder1,folder2,folder,fid,WarningMessage] = clean_mkr_part1(markers1,markers2,folder1,folder2,folder,fid,WarningMessage,DINtype)

for i=1:2
    switch i
        case 1; mkr_file = markers1;
        case 2; mkr_file = markers2;
    end
    mkr = dual_read_mkr(mkr_file,DINtype);
    
    %% Non-numerical
    if any(isnan(mkr(:,2)))
        WarningMessage(i,1)=1;
        warning('Discarding non-numerical markers, check visually in the file %s', mkr_file)
        mkr(isnan(mkr(:,2)),:)=[];
    end
    
    NonNumMark{i} = mkr((isnan(mkr(:,2))),1);
    
    %% Multiple labels "1"
    multi_ones = find(mkr(:,2)==1);
%     multi_ones(multi_ones==1) = []; % spare first one
    if ~isempty(multi_ones)
        multi_ones(1) = []; % spare first one #RM, 2017-08-17
    end
    if ~isempty(multi_ones)
        % should not be done before unwrapping! (reloading markers for
        % unwrapping will be a workaround for this)
        mkr(multi_ones,:)=[];
        WarningMessage(i,2)=1;
        warning('Discarding marker(s) that appeared as "1" at rows: ')
        disp(multi_ones)
    end
    % if labels were unwrapped, it satisfies the first condition and it
    % will overwrite the incorrectly "multi-ones" detected before unwrapping
    % and if labels did not need to be unwrapped, the second condition
    % is satisfied and will not be overwritten later because Task will
    % never be 'donotreload' without being 'unwrap' before #RM
    MultiOnes{i} = multi_ones;
    
    %% Label doublets and gaps
%     delta = [0 diff(mkr(:,2))]; %#RM-2017-08-14
    delta = diff(mkr(:,2));
    % !! What follows has to be done before unwrapping, otherwise unwrapping
    % will later generated gaps in the markers (when delta>1) !!
    % (RM,2017-08-16)
%     if ~strcmpi(Task,'unwrap')
        
%         if sum(delta==0) > 1 %#RM-2017-08-14
        if sum(delta==0) > 0
            % in this case, it is a "single neighbouring doublet" case
            mkr(delta==0,:)=[]; % get rid of first of each doublet
            WarningMessage(i,3)=1;
            warning('Discarding %d marker doublets',sum(delta==0))
        end
        
        DoubletsGaps{i} = delta;
        
        delta = diff(mkr(:,2));
        % separated delta==0 and delta>1 because both delta == 0 and delta > 1 can be true
        if any(delta>1) % but if delta < 1 then it can only means that unwrapping is not done or that there are doublets (and the latter are removed by detecting delta > 1 concurrently with delta < 1
            
            %#RM-2017-08-16: Could be a "shifted doublet"!!
            SingleShiftedDoubletTmp = find(delta>1);
            % If after the suspicious marker label diff gets negative and
            % closer to zero by 1, then it is a "single shifted doublet",
            % and we remove it:
            % ... But before we check that these these are not
            SingleShiftedDoublet = nan(size(SingleShiftedDoubletTmp));
            for Nssd = 1:length(SingleShiftedDoubletTmp)
                if SingleShiftedDoubletTmp(Nssd)==length(delta)
                    % - Reached the end?
                    % - Simply remove it, useless to get a single marker
                    % after a gap...
                    SingleShiftedDoublet(Nssd)=length(delta)+1;
                else
                    % if diff becomes negative after the gap, it was a
                    % "single shifted doublet"
%                     if -(delta(SingleShiftedDoubletTmp(Nssd))-1)==delta(SingleShiftedDoubletTmp(Nssd)+1)
%                         SingleShiftedDoublet(Nssd)=SingleShiftedDoubletTmp(Nssd);
                    if delta(SingleShiftedDoubletTmp(Nssd)+1)<0
                        SingleShiftedDoublet(Nssd)=SingleShiftedDoubletTmp(Nssd)+1;
                    else
                        SingleShiftedDoublet(Nssd)=nan;
                    end
                end
            end
            SingleShiftedDoublet = SingleShiftedDoublet(~isnan(SingleShiftedDoublet));
            
            if size(SingleShiftedDoublet,1)>=size(mkr,1)/4
                % ...too much to remove...
                % maybe the serie of markers contains an extra set of
                % marker, and it might start with the wrong one
                
                % ONLY DO THIS TRICKY PART IF THE SITUATION IS AMBIGUOUS !!!!
                if i==1 % OK to load second file, anyhow it will be reloaded from 0 later
                    other_mkr = dual_read_mkr(markers2,DINtype);
                elseif i==2
                    % first file already processed, do not reload it from 0
                    % to compare with second file ! %#RM, 2017-08-17
                    other_mkr = mkr1;
                end
                
                ExtraIdx = match_vectors(setdiff(mkr(:,2),other_mkr(:,2)),mkr(:,2),1);
                if isempty(ExtraIdx)
                    [~,BigMkr] = find([size(mkr,1),size(other_mkr,1)]==max([size(mkr,1),size(other_mkr,1)]));
                    if length(BigMkr)<2
                            TempWarpedDoublets1 = find(cellfun(@length,match_vectors(other_mkr(:,2),mkr(:,2),1))>1);
                            TempWarpedDoublets2 = match_vectors(other_mkr(:,2),mkr(:,2),1);
                            DeleteThisWarpedDoublet = nan(size(TempWarpedDoublets1));
                            for twd = 1:length(TempWarpedDoublets1)
                                DeleteThisWarpedDoublet(twd) = TempWarpedDoublets2{TempWarpedDoublets1(twd)}(2);
                            end
                        mkr(DeleteThisWarpedDoublet,:)=[];
                    end
                elseif isnumeric(ExtraIdx)
                    mkr(ExtraIdx,:)=[];
                elseif iscell(ExtraIdx)
                    mkr(cell2mat(ExtraIdx),:)=[];
                end
                WarningMessage(i,5)=1;
                warning('Discarding %d extra markers',sum(delta<1))
                
            elseif size(SingleShiftedDoublet,1)==0
                % ...nothing to remove...
                % which means we are facing "true" gaps
                WarningMessage(i,4)=1;
                warning('Gaps in the markers, need to check visually file %s', mkr_file)
                disp(mkr(delta>1,2))
            else
                WarningMessage(i,5)=1;
                warning('Discarding %d extra markers',sum(delta<1))
                % Remove the "single shifted doublets"
                mkr(SingleShiftedDoublet,:)=[];
            end
            
        end
%     end
    
    ExtraMarkersDiscarded{i} = delta;

    switch i
        case 1; mkr1 = mkr;
        case 2; mkr2 = mkr;
    end
end

%% Markers missing in one of the two files
%#RM, 2017-08-16, NB: markers at beginning and end of files can be continuous
% and not flagged as suspicious by previous checks but they still miss in
% the other file because recordings do not stop exactly at the same time

[TempIntruders,idx] = setdiff(mkr1(:,2),mkr2(:,2));
if ~isempty(idx)
    % RM, Aug 2017 line below works except if there are multiple occurences
    % of these markers present only in 1 out of 2 files! (see help setdiff.m)
    %     mkr1(idx,:)=[];
    idx = 0;
    for mm = 1:length(TempIntruders)
        mkr1(mkr1(:,2)==TempIntruders(mm),:)=[];
        idx = idx+sum(mkr1(:,2)==TempIntruders(mm));
    end
    WarningMessage(1,6)=1;
    warning('Discarding %d markers from file 1 that were not present in file 2',numel(idx))
end
icEEGbutNOThdEEG = idx;
[TempIntruders,idx] = setdiff(mkr2(:,2),mkr1(:,2));
if ~isempty(idx)
    % RM, Aug 2017 line below works except if there are multiple occurences
    % of these markers present only in 1 out of 2 files! (see help setdiff.m)
    %     mkr2(idx,:)=[];
    idx = 0;
    for mm = 1:length(TempIntruders)
        mkr2(mkr2(:,2)==TempIntruders(mm),:)=[];
        idx = idx+sum(mkr2(:,2)==TempIntruders(mm));
    end
    WarningMessage(2,7)=1;
    warning('Discarding %d markers from file 2 that were not present in file 1',numel(idx))
end
hdEEGbutNOTicEEG = idx;

end

function [mkr, delta, mkr1, mkr2, DoubletsGaps, ExtraMarkersDiscarded, markers1,markers2,folder1,folder2,folder,fid,WarningMessage] = clean_mkr_unwrap(markers1,markers2,folder1,folder2,folder,fid,WarningMessage,DINtype)

for i=1:2
    switch i
        case 1; mkr_file = markers1;
        case 2; mkr_file = markers2;
    end
    mkr = dual_read_mkr(mkr_file,DINtype);
    
    %% Label doublets and gaps
%     delta = [0 diff(mkr(:,2))]; %#RM-2017-08-14
    delta = diff(mkr(:,2));
    % !! What follows has to be done before unwrapping, otherwise unwrapping
    % will later generated gaps in the markers (when delta>1) !!
    % (RM,2017-08-16)
%     if ~strcmpi(Task,'unwrap')
        
%         if sum(delta==0) > 1 %#RM-2017-08-14
        if sum(delta==0) > 0
            % in this case, it is a "single neighbouring doublet" case
            mkr(delta==0,:)=[]; % get rid of first of each doublet
            WarningMessage(i,3)=1;
            warning('Discarding %d marker doublets',sum(delta==0))
        end
        
        DoubletsGaps{i} = delta;
        
        delta = diff(mkr(:,2));
        % separated delta==0 and delta>1 because both delta == 0 and delta > 1 can be true
        if any(delta>1) % but if delta < 1 then it can only means that unwrapping is not done or that there are doublets (and the latter are removed by detecting delta > 1 concurrently with delta < 1
            
            %#RM-2017-08-16: Could be a "shifted doublet"!!
            SingleShiftedDoubletTmp = find(delta>1);
            % If after the suspicious marker label diff gets negative and
            % closer to zero by 1, then it is a "single shifted doublet",
            % and we remove it:
            % ... But before we check that these these are not
            SingleShiftedDoublet = nan(size(SingleShiftedDoubletTmp));
            for Nssd = 1:length(SingleShiftedDoubletTmp)
                if SingleShiftedDoubletTmp(Nssd)==length(delta)
                    % - Reached the end?
                    % - Simply remove it, useless to get a single marker
                    % after a gap...
                    SingleShiftedDoublet(Nssd)=length(delta)+1;
                else
                    % if diff becomes negative after the gap, it was a
                    % "single shifted doublet"
%                     if -(delta(SingleShiftedDoubletTmp(Nssd))-1)==delta(SingleShiftedDoubletTmp(Nssd)+1)
%                         SingleShiftedDoublet(Nssd)=SingleShiftedDoubletTmp(Nssd);
                    if delta(SingleShiftedDoubletTmp(Nssd)+1)<0
                        SingleShiftedDoublet(Nssd)=SingleShiftedDoubletTmp(Nssd)+1;
                    else
                        SingleShiftedDoublet(Nssd)=nan;
                    end
                end
            end
            SingleShiftedDoublet = SingleShiftedDoublet(~isnan(SingleShiftedDoublet));
            
            if size(SingleShiftedDoublet,1)>=size(mkr,1)/4
                % ...too much to remove...
                % maybe the serie of markers contains an extra set of
                % marker, and it might start with the wrong one
                
                % ONLY DO THIS TRICKY PART IF THE SITUATION IS AMBIGUOUS !!!!
                if i==1 % OK to load second file, anyhow it will be reloaded from 0 later
                    other_mkr = dual_read_mkr(markers2,DINtype);
                elseif i==2
                    % first file already processed, do not reload it from 0
                    % to compare with second file ! %#RM, 2017-08-17
                    other_mkr = mkr1;
                end
                
                ExtraIdx = match_vectors(setdiff(mkr(:,2),other_mkr(:,2)),mkr(:,2),1);
                if isempty(ExtraIdx)
                    [~,BigMkr] = find([size(mkr,1),size(other_mkr,1)]==max([size(mkr,1),size(other_mkr,1)]));
                    if length(BigMkr)<2
                            TempWarpedDoublets1 = find(cellfun(@length,match_vectors(other_mkr(:,2),mkr(:,2),1))>1);
                            TempWarpedDoublets2 = match_vectors(other_mkr(:,2),mkr(:,2),1);
                            DeleteThisWarpedDoublet = nan(size(TempWarpedDoublets1));
                            for twd = 1:length(TempWarpedDoublets1)
                                DeleteThisWarpedDoublet(twd) = TempWarpedDoublets2{TempWarpedDoublets1(twd)}(2);
                            end
                        mkr(DeleteThisWarpedDoublet,:)=[];
                    end
                elseif isnumeric(ExtraIdx)
                    mkr(ExtraIdx,:)=[];
                elseif iscell(ExtraIdx)
                    mkr(cell2mat(ExtraIdx),:)=[];
                end
                WarningMessage(i,5)=1;
                warning('Discarding %d extra markers',sum(delta<1))
                
            elseif size(SingleShiftedDoublet,1)==0
                % ...nothing to remove...
                % which means we are facing "true" gaps
                WarningMessage(i,4)=1;
                warning('Gaps in the markers, need to check visually file %s', mkr_file)
                disp(mkr(delta>1,2))
            else
                WarningMessage(i,5)=1;
                warning('Discarding %d extra markers',sum(delta<1))
                % Remove the "single shifted doublets"
                mkr(SingleShiftedDoublet,:)=[];
            end
            
        end
%     end
    
    ExtraMarkersDiscarded{i} = delta;
    
    %% /////// LABEL UNWRAPPING ///////
    if any(diff(mkr(:,2))<0)
        WarningMessage(i,10)=1;
    end
    IdxSteps = find(diff(mkr(:,2))<0);
    mkrBKP = mkr;
    for nn = 1:length(IdxSteps)
        mkr(IdxSteps(nn)+1:end,2) = mkr(IdxSteps(nn)+1:end,2) + mkrBKP(IdxSteps(nn),2);
    end
    
    switch i
        case 1; mkr1 = mkr;
        case 2; mkr2 = mkr;
    end
end

end

function [mkr, multi_ones, MultiOnes, icEEGbutNOThdEEG, hdEEGbutNOTicEEG, mkr1, mkr2, idx, markers1,markers2,folder1,folder2,folder,fid,WarningMessage] = clean_mkr_part2(markers1,markers2,folder1,folder2,folder,fid,WarningMessage,varargin)

for i=1:2
    
    switch i
        case 1
            mkr = varargin{1}; % i.e. mkr1temp
            mkr_file = markers1;
        case 2
            mkr = varargin{2}; % i.e. mkr2temp
            mkr_file = markers2;
    end
    
    %% Multiple labels "1"
    multi_ones = find(mkr(:,2)==1);
%     multi_ones(multi_ones==1) = []; % spare first one
    if ~isempty(multi_ones)
        multi_ones(1) = []; % spare first one #RM, 2017-08-17
    end
    if ~isempty(multi_ones)
        mkr(multi_ones,:)=[];
        WarningMessage(i,2)=1;
        warning('Discarding marker(s) that appeared as "1" at rows: ')
        disp(multi_ones)
    else
        WarningMessage(i,2)=0; % see comment below...
    end
    
    % if labels were unwrapped, it satisfies the first condition and it
    % will overwrite the incorrectly "multi-ones" detected before unwrapping
    % and if labels did not need to be unwrapped, the second condition
    % is satisfied and will not be overwritten later because Task will
    % never be 'donotreload' without being 'unwrap' before #RM
    MultiOnes{i} = multi_ones;
    
    %% LABEL DOUBLETS - LIGHT VERSION
    delta = diff(mkr(:,2));
    if sum(delta==0) > 0
        % in this case, it is a "single neighbouring doublet" case
        mkr(delta==0,:)=[]; % get rid of first of each doublet
        WarningMessage(i,3)=1;
        warning('Discarding %d marker doublets',sum(delta==0))
    end
    DoubletsGaps{i} = delta;
    ExtraMarkersDiscarded{i} = delta;
    
    switch i
        case 1; mkr1 = mkr;
        case 2; mkr2 = mkr;
    end
end

%% Markers missing in one of the two files
%#RM, 2017-08-16, NB: markers at beginning and end of files can be continuous
% and not flagged as suspicious by previous checks but they still miss in
% the other file because recordings do not stop exactly at the same time

[TempIntruders,idx] = setdiff(mkr1(:,2),mkr2(:,2));
if ~isempty(idx)
    % RM, Aug 2017 line below works except if there are multiple occurences
    % of these markers present only in 1 out of 2 files! (see help setdiff.m)
    %     mkr1(idx,:)=[];
    idx = 0;
    for mm = 1:length(TempIntruders)
        mkr1(mkr1(:,2)==TempIntruders(mm),:)=[];
        idx = idx+sum(mkr1(:,2)==TempIntruders(mm));
    end
    WarningMessage(1,6)=1;
    warning('Discarding %d markers from file 1 that were not present in file 2',numel(idx))
else
    WarningMessage(1,6)=0;
end
icEEGbutNOThdEEG = idx;
[TempIntruders,idx] = setdiff(mkr2(:,2),mkr1(:,2));
if ~isempty(idx)
    % RM, Aug 2017 line below works except if there are multiple occurences
    % of these markers present only in 1 out of 2 files! (see help setdiff.m)
    %     mkr2(idx,:)=[];
    idx = 0;
    for mm = 1:length(TempIntruders)
        mkr2(mkr2(:,2)==TempIntruders(mm),:)=[];
        idx = idx+sum(mkr2(:,2)==TempIntruders(mm));
    end
    WarningMessage(2,7)=1;
    warning('Discarding %d markers from file 2 that were not present in file 1',numel(idx))
else
    WarningMessage(2,7)=0;
end
hdEEGbutNOTicEEG = idx;

end

function [mkr, NonNumMark, multi_ones, MultiOnes, delta, DoubletsGaps, ExtraMarkersDiscarded, icEEGbutNOThdEEG, hdEEGbutNOTicEEG, mkr1, mkr2, idx, markers1,markers2,folder1,folder2,folder,fid,WarningMessage] = clean_mkr_part1_but_deunwrap_before(markers1,markers2,folder1,folder2,folder,fid,WarningMessage,Period,DINtype)

for i=1:2
    switch i
        case 1; mkr_file = markers1;
        case 2; mkr_file = markers2;
    end
    mkr = dual_read_mkr(mkr_file,DINtype);
    
    %% Non-numerical
    if any(isnan(mkr(:,2)))
        WarningMessage(i,1)=1;
        warning('Discarding non-numerical markers, check visually in the file %s', mkr_file)
        mkr(isnan(mkr(:,2)),:)=[];
    end
    
    NonNumMark{i} = mkr((isnan(mkr(:,2))),1);
    
    %% Multiple labels "1"
    multi_ones = find(mkr(:,2)==1);
%     multi_ones(multi_ones==1) = []; % spare first one
    if ~isempty(multi_ones)
        multi_ones(1) = []; % spare first one #RM, 2017-08-17
    end
    if ~isempty(multi_ones)
        % should not be done before unwrapping! (reloading markers for
        % unwrapping will be a workaround for this)
        mkr(multi_ones,:)=[];
        WarningMessage(i,2)=1;
        warning('Discarding marker(s) that appeared as "1" at rows: ')
        disp(multi_ones)
    end
    % if labels were unwrapped, it satisfies the first condition and it
    % will overwrite the incorrectly "multi-ones" detected before unwrapping
    % and if labels did not need to be unwrapped, the second condition
    % is satisfied and will not be overwritten later because Task will
    % never be 'donotreload' without being 'unwrap' before #RM
    MultiOnes{i} = multi_ones;
    
    %% Label doublets and gaps
%     delta = [0 diff(mkr(:,2))]; %#RM-2017-08-14
    delta = diff(mkr(:,2));
    % !! What follows has to be done before unwrapping, otherwise unwrapping
    % will later generated gaps in the markers (when delta>1) !!
    % (RM,2017-08-16)
%     if ~strcmpi(Task,'unwrap')
        
%         if sum(delta==0) > 1 %#RM-2017-08-14
        if sum(delta==0) > 0
            % in this case, it is a "single neighbouring doublet" case
            mkr(delta==0,:)=[]; % get rid of first of each doublet
            WarningMessage(i,3)=1;
            warning('Discarding %d marker doublets',sum(delta==0))
        end
        
        DoubletsGaps{i} = delta;
        
        delta = diff(mkr(:,2));
        % separated delta==0 and delta>1 because both delta == 0 and delta > 1 can be true
        if any(delta>1) % but if delta < 1 then it can only means that unwrapping is not done or that there are doublets (and the latter are removed by detecting delta > 1 concurrently with delta < 1
            
            %#RM-2017-08-16: Could be a "shifted doublet"!!
            SingleShiftedDoubletTmp = find(delta>1);
            % If after the suspicious marker label diff gets negative and
            % closer to zero by 1, then it is a "single shifted doublet",
            % and we remove it:
            % ... But before we check that these these are not
            SingleShiftedDoublet = nan(size(SingleShiftedDoubletTmp));
            for Nssd = 1:length(SingleShiftedDoubletTmp)
                if SingleShiftedDoubletTmp(Nssd)==length(delta)
                    % - Reached the end?
                    % - Simply remove it, useless to get a single marker
                    % after a gap...
                    SingleShiftedDoublet(Nssd)=length(delta)+1;
                else
                    % if diff becomes negative after the gap, it was a
                    % "single shifted doublet"
%                     if -(delta(SingleShiftedDoubletTmp(Nssd))-1)==delta(SingleShiftedDoubletTmp(Nssd)+1)
%                         SingleShiftedDoublet(Nssd)=SingleShiftedDoubletTmp(Nssd);
                    if delta(SingleShiftedDoubletTmp(Nssd)+1)<0
                        SingleShiftedDoublet(Nssd)=SingleShiftedDoubletTmp(Nssd)+1;
                    else
                        SingleShiftedDoublet(Nssd)=nan;
                    end
                end
            end
            SingleShiftedDoublet = SingleShiftedDoublet(~isnan(SingleShiftedDoublet));
            
            if size(SingleShiftedDoublet,1)>=size(mkr,1)/4
                % ...too much to remove...
                % maybe the serie of markers contains an extra set of
                % marker, and it might start with the wrong one
                
                % ONLY DO THIS TRICKY PART IF THE SITUATION IS AMBIGUOUS !!!!
                if i==1 % OK to load second file, anyhow it will be reloaded from 0 later
                    other_mkr = dual_read_mkr(markers2,DINtype);
                elseif i==2
                    % first file already processed, do not reload it from 0
                    % to compare with second file ! %#RM, 2017-08-17
                    other_mkr = mkr1;
                end
                
                ExtraIdx = match_vectors(setdiff(mkr(:,2),other_mkr(:,2)),mkr(:,2),1);
                if isempty(ExtraIdx)
                    [~,BigMkr] = find([size(mkr,1),size(other_mkr,1)]==max([size(mkr,1),size(other_mkr,1)]));
                    if length(BigMkr)<2
                            TempWarpedDoublets1 = find(cellfun(@length,match_vectors(other_mkr(:,2),mkr(:,2),1))>1);
                            TempWarpedDoublets2 = match_vectors(other_mkr(:,2),mkr(:,2),1);
                            DeleteThisWarpedDoublet = nan(size(TempWarpedDoublets1));
                            for twd = 1:length(TempWarpedDoublets1)
                                DeleteThisWarpedDoublet(twd) = TempWarpedDoublets2{TempWarpedDoublets1(twd)}(2);
                            end
                        mkr(DeleteThisWarpedDoublet,:)=[];
                    end
                elseif isnumeric(ExtraIdx)
                    mkr(ExtraIdx,:)=[];
                elseif iscell(ExtraIdx)
                    mkr(cell2mat(ExtraIdx),:)=[];
                end
                WarningMessage(i,5)=1;
                warning('Discarding %d extra markers',sum(delta<1))
                
            elseif size(SingleShiftedDoublet,1)==0
                % ...nothing to remove...
                % which means we are facing "true" gaps
                WarningMessage(i,4)=1;
                warning('Gaps in the markers, need to check visually file %s', mkr_file)
                disp(mkr(delta>1,2))
            else
                WarningMessage(i,5)=1;
                warning('Discarding %d extra markers',sum(delta<1))
                % Remove the "single shifted doublets"
                mkr(SingleShiftedDoublet,:)=[];
            end
            
        end
%     end
    
    ExtraMarkersDiscarded{i} = delta;

    switch i
        case 1; mkr1 = mkr;
        case 2; mkr2 = mkr;
    end
end

% De-unwrap if needed...
if (any(diff(mkr1(:,2))<0) && all(diff(mkr2(:,2))>=0)) || (any(diff(mkr2(:,2))<0) && all(diff(mkr1(:,2))>=0))
    
    if (any(diff(mkr1(:,2))<0) && all(diff(mkr2(:,2))>=0))
        
        for nn = 1:length(mkr2(:,1))
            if mod(mkr2(nn,2),Period)~=0
                mkr2(nn,2)=mod(mkr2(nn,2),Period);
            else
                mkr2(nn,2)=Period; % WORKS BECAUSE THERE IS ONLY 1 PERIOD!!
            end
        end
        
    elseif (any(diff(mkr2(:,2))<0) && all(diff(mkr1(:,2))>=0))
        
        for nn = 1:length(mkr1(:,1))
            if mod(mkr1(nn,2),Period)~=0
                mkr1(nn,2)=mod(mkr1(nn,2),Period);
            else
                mkr1(nn,2)=Period; % WORKS BECAUSE THERE IS ONLY 1 PERIOD!!
            end
        end
        
    end
    
end

%% Markers missing in one of the two files
%#RM, 2017-08-16, NB: markers at beginning and end of files can be continuous
% and not flagged as suspicious by previous checks but they still miss in
% the other file because recordings do not stop exactly at the same time

[TempIntruders,idx] = setdiff(mkr1(:,2),mkr2(:,2));
if ~isempty(idx)
    % RM, Aug 2017 line below works except if there are multiple occurences
    % of these markers present only in 1 out of 2 files! (see help setdiff.m)
    %     mkr1(idx,:)=[];
    idx = 0;
    for mm = 1:length(TempIntruders)
        mkr1(mkr1(:,2)==TempIntruders(mm),:)=[];
        idx = idx+sum(mkr1(:,2)==TempIntruders(mm));
    end
    WarningMessage(1,6)=1;
    warning('Discarding %d markers from file 1 that were not present in file 2',numel(idx))
end
icEEGbutNOThdEEG = idx;
[TempIntruders,idx] = setdiff(mkr2(:,2),mkr1(:,2));
if ~isempty(idx)
    % RM, Aug 2017 line below works except if there are multiple occurences
    % of these markers present only in 1 out of 2 files! (see help setdiff.m)
    %     mkr2(idx,:)=[];
    idx = 0;
    for mm = 1:length(TempIntruders)
        mkr2(mkr2(:,2)==TempIntruders(mm),:)=[];
        idx = idx+sum(mkr2(:,2)==TempIntruders(mm));
    end
    WarningMessage(2,7)=1;
    warning('Discarding %d markers from file 2 that were not present in file 1',numel(idx))
end
hdEEGbutNOTicEEG = idx;

end

function [mkr, delta, mkr1, mkr2, DoubletsGaps, ExtraMarkersDiscarded, markers1,markers2,folder1,folder2,folder,fid,WarningMessage] = clean_mkr_deunwrap_and_unwrap(markers1,markers2,folder1,folder2,folder,fid,WarningMessage,Period,DINtype)

for i=1:2
    switch i
        case 1; mkr_file = markers1;
        case 2; mkr_file = markers2;
    end
    mkr = dual_read_mkr(mkr_file,DINtype);
    
    %% Label doublets and gaps
%     delta = [0 diff(mkr(:,2))]; %#RM-2017-08-14
    delta = diff(mkr(:,2));
    % !! What follows has to be done before unwrapping, otherwise unwrapping
    % will later generated gaps in the markers (when delta>1) !!
    % (RM,2017-08-16)
%     if ~strcmpi(Task,'unwrap')
        
%         if sum(delta==0) > 1 %#RM-2017-08-14
        if sum(delta==0) > 0
            % in this case, it is a "single neighbouring doublet" case
            mkr(delta==0,:)=[]; % get rid of first of each doublet
            WarningMessage(i,3)=1;
            warning('Discarding %d marker doublets',sum(delta==0))
        end
        
        DoubletsGaps{i} = delta;
        
        delta = diff(mkr(:,2));
        % separated delta==0 and delta>1 because both delta == 0 and delta > 1 can be true
        if any(delta>1) % but if delta < 1 then it can only means that unwrapping is not done or that there are doublets (and the latter are removed by detecting delta > 1 concurrently with delta < 1
            
            %#RM-2017-08-16: Could be a "shifted doublet"!!
            SingleShiftedDoubletTmp = find(delta>1);
            % If after the suspicious marker label diff gets negative and
            % closer to zero by 1, then it is a "single shifted doublet",
            % and we remove it:
            % ... But before we check that these these are not
            SingleShiftedDoublet = nan(size(SingleShiftedDoubletTmp));
            for Nssd = 1:length(SingleShiftedDoubletTmp)
                if SingleShiftedDoubletTmp(Nssd)==length(delta)
                    % - Reached the end?
                    % - Simply remove it, useless to get a single marker
                    % after a gap...
                    SingleShiftedDoublet(Nssd)=length(delta)+1;
                else
                    % if diff becomes negative after the gap, it was a
                    % "single shifted doublet"
%                     if -(delta(SingleShiftedDoubletTmp(Nssd))-1)==delta(SingleShiftedDoubletTmp(Nssd)+1)
%                         SingleShiftedDoublet(Nssd)=SingleShiftedDoubletTmp(Nssd);
                    if delta(SingleShiftedDoubletTmp(Nssd)+1)<0
                        SingleShiftedDoublet(Nssd)=SingleShiftedDoubletTmp(Nssd)+1;
                    else
                        SingleShiftedDoublet(Nssd)=nan;
                    end
                end
            end
            SingleShiftedDoublet = SingleShiftedDoublet(~isnan(SingleShiftedDoublet));
            
            if size(SingleShiftedDoublet,1)>=size(mkr,1)/4
                % ...too much to remove...
                % maybe the serie of markers contains an extra set of
                % marker, and it might start with the wrong one
                
                % ONLY DO THIS TRICKY PART IF THE SITUATION IS AMBIGUOUS !!!!
                if i==1 % OK to load second file, anyhow it will be reloaded from 0 later
                    other_mkr = dual_read_mkr(markers2,DINtype);
                elseif i==2
                    % first file already processed, do not reload it from 0
                    % to compare with second file ! %#RM, 2017-08-17
                    other_mkr = mkr1;
                end
                
                ExtraIdx = match_vectors(setdiff(mkr(:,2),other_mkr(:,2)),mkr(:,2),1);
                if isempty(ExtraIdx)
                    [~,BigMkr] = find([size(mkr,1),size(other_mkr,1)]==max([size(mkr,1),size(other_mkr,1)]));
                    if length(BigMkr)<2
                            TempWarpedDoublets1 = find(cellfun(@length,match_vectors(other_mkr(:,2),mkr(:,2),1))>1);
                            TempWarpedDoublets2 = match_vectors(other_mkr(:,2),mkr(:,2),1);
                            DeleteThisWarpedDoublet = nan(size(TempWarpedDoublets1));
                            for twd = 1:length(TempWarpedDoublets1)
                                DeleteThisWarpedDoublet(twd) = TempWarpedDoublets2{TempWarpedDoublets1(twd)}(2);
                            end
                        mkr(DeleteThisWarpedDoublet,:)=[];
                    end
                elseif isnumeric(ExtraIdx)
                    mkr(ExtraIdx,:)=[];
                elseif iscell(ExtraIdx)
                    mkr(cell2mat(ExtraIdx),:)=[];
                end
                WarningMessage(i,5)=1;
                warning('Discarding %d extra markers',sum(delta<1))
                
            elseif size(SingleShiftedDoublet,1)==0
                % ...nothing to remove...
                % which means we are facing "true" gaps
                WarningMessage(i,4)=1;
                warning('Gaps in the markers, need to check visually file %s', mkr_file)
                disp(mkr(delta>1,2))
            else
                WarningMessage(i,5)=1;
                warning('Discarding %d extra markers',sum(delta<1))
                % Remove the "single shifted doublets"
                mkr(SingleShiftedDoublet,:)=[];
            end
            
        end
%     end
    
    ExtraMarkersDiscarded{i} = delta;
    
    %% /////// LABEL DE-UNWRAPPING AND RE-UNWRAPPING ///////
    
    % De-unwrap if needed...
    for nn = 1:length(mkr(:,1))
        if mod(mkr(nn,2),Period)~=0
            mkr(nn,2)=mod(mkr(nn,2),Period);
        else
            mkr(nn,2)=Period; % WORKS BECAUSE THERE IS ONLY 1 PERIOD!!
        end
    end
    
    %% /////// LABEL UNWRAPPING ///////
    if any(diff(mkr(:,2))<0)
        WarningMessage(i,10)=1;
    end
    IdxSteps = find(diff(mkr(:,2))<0);
    mkrBKP = mkr;
    for nn = 1:length(IdxSteps)
        mkr(IdxSteps(nn)+1:end,2) = mkr(IdxSteps(nn)+1:end,2) + mkrBKP(IdxSteps(nn),2);
    end
    
    switch i
        case 1; mkr1 = mkr;
        case 2; mkr2 = mkr;
    end
end

end
