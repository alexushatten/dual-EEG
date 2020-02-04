function spike_aligner_wrapper( MrkFile, MarkersOfInterest, Fs )
% spike_aligner_wrapper: wrapper to align individual spikes with average
% and write new .mrk file
%
% spike_aligner_wrapper( MrkFile, MarkersOfInterest, Fs )
%
%  Inputs
% --------
% MrkFile: 
% MarkersOfInterest: 
% Fs: sampling rate
%
%  Outputs
% ---------
% A backup of original .mrk file with suffix "_BKP_before_align"
%
% A modified .mrk file with shifted time stamps for spikes (other markers
% are preserved)
%
% A .mat file with original epochs and adjustments applied
%
%-------------------------------------------------------------------------
% NB:
%
% /!\ TO EXCLUDE SPIKES FROM AVERAGE (append "_excluded" to the marker
% label), MARK THEM AS "TO CHECK" DURING INTERACTIVE REVIEW PROCESS!
%
% /!\ DO NOT FORGET TO (BKP AND) RE-RUN AVERAGING AFTER EXECUTING THIS FUNCTION!
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, September 2018
%-------------------------------------------------------------------------

% here we load results for bipolar montage
% (we do not use monopolar montage with physical reference as spike marking
% was done using bipolar montage!)
for s = 1:length(MarkersOfInterest)
    MarkerType = MarkersOfInterest{s};
    ResultsFiles{s} = [pwd,filesep,'spk_avg_icEEG_',regexprep(regexprep(regexprep(MarkerType,'"',''),'0',''),'_',''),'_Bipolar.mat'];
    if exist(ResultsFiles{s},'file')~=2
        error(['File ',ResultsFiles{s},' was not found on disk!'])
    end
end; ResultsFiles = ResultsFiles';

% this has to be done outside of loop for each spike type, otherwise the
% original version will get lost when there are multiple spike types!
for f = 1:length(MrkFile)
    copyfile(MrkFile{f},spm_file(MrkFile{f},'suffix','_BKP_before_align'));
end

% review and adjust each spike type
for s = 1:length(MarkersOfInterest)
    load(ResultsFiles{s},'BCEEGic','icLabels','SpikeTimings')
    
    % interactive review
    [Shifts, ToCheck] = spike_aligner(BCEEGic,Fs,icLabels);
    
    % calculate new time stamps & find spikes to exclude from average
    NewSpikeTimings = SpikeTimings;
    for f = 1:size(SpikeTimings(:,s),1)
    Count = 0;
        for sp = 1:length(SpikeTimings{f,s})
            Count = Count + 1;
            if ToCheck(Count)
                NewSpikeTimings{f,s}(Count) = nan;
            else
                NewSpikeTimings{f,s}(Count) = NewSpikeTimings{f,s}(Count)+Shifts(Count);
            end
        end
    end
    
    % write new time stamps
    for f = 1:size(SpikeTimings(:,s),1)
        [MarkerTimeStart, MarkerTimeEnd, MarkerLabel] = read_mrk_file(MrkFile{f});
        NewMarkerTimeStart = MarkerTimeStart;
        NewMarkerTimeEnd = MarkerTimeEnd;
        NewMarkerLabel = MarkerLabel;
        TempMatch = match_vectors(SpikeTimings{f,s},MarkerTimeStart,0);
        if isnumeric(TempMatch) % just to be safe
            NewMarkerTimeStart(TempMatch) = NewSpikeTimings{f,s};
            NewMarkerTimeEnd(TempMatch) = NewSpikeTimings{f,s};
        end
        for tf = 1:length(NewMarkerTimeStart)
            if isnan(NewMarkerTimeStart(tf))
                NewMarkerLabel{tf} = ['"',regexprep(NewMarkerLabel{tf},'"',''),'_excluded"'];
                NewMarkerTimeStart(tf) = MarkerTimeStart(tf);
                NewMarkerTimeEnd(tf) = MarkerTimeEnd(tf);
            end
        end
        write_mrk_file_Cartool(MrkFile{f}, NewMarkerTimeStart, NewMarkerTimeEnd, NewMarkerLabel);
    end
    
    % save shifts
    save([pwd,filesep,'spike_aligner_adjustments_',regexprep(regexprep(regexprep(MarkerType,'"',''),'0',''),'_',''),'.mat'],'Shifts','ToCheck','BCEEGic','icLabels','SpikeTimings','NewSpikeTimings');
end

end

