function [Bad_Channels, Good_Channels, Bad_Segments, Good_Segments, Metrics, Results, Settings] = rm_artefacts(EEG, SamplingFreq, RefChan, MontageInfo)
% rm_artefacts: rm(ARTEFACTS) is an artefacts detection pipeline,
% implementing recursive stripped down version of FASTER algorithm,
% combined with GFP and Dissimilarity measures, as well as a my new metric
% that targets the detection of single (or few) electrodes with outlier
% signal, such as e.g. electrode pops. Regarding the FASTER part, compared
% to the original methods, the z-threshold values were raised up to avoid
% rejecting good EEG segments such as the ones containing a lot of alpha
% (we do not want to excise portions of EEG containing strong alpha waves,
% if this has to be removed, this is more something to be done at the ICA
% stage), and uses only the 1st and 2nd stages of the original pipeline
% (bad channels and bad segments detection, across the whole recording, and
% for each epoch) with some modifications. Compared to the original method,
% we would prefer to feed the ICA with as-much-as posible continuous EEG
% traces (although this can be debated from a methodological perspective,
% this is mostly for practical reasons, e.g. avoid having to perform again
% the ICA if the epochs change), therefore this pipeline generates
% overlapping epochs over time and rejects time points which are
% consistenly detected as bad in the segments dimension.
% rm(ARTEFACTS) does not implement the ICA stage of FASTER, and physiological
% artefacts such as eye blinks, eye movements and cardiac-related artefacts,
% are not detected by rm(ARTEFACTS) and shall be removed at later ICA stage.
% The purpose of rm(ARTEFACTS) is therefore to give the ICA a "clean" EEG
% signal without major non-physiological artefacts. It has been tested
% on 257-channels scalp EEG data (EGI).
%
% [Bad_Channels, Good_Channels, ...
%  Bad_Segments, Good_Segments, ...
%  Metrics, Results, Settings] = rm_artefacts(EEG, SamplingFreq, ...
%                                 RefChan, MontageInfo)
%
%  Inputs
% --------
% EEG: [channel x time] EEG traces, ALREADY FILTERED!
% SamplingFreq: [1 x 1] integer, sampling frequency
% RefChan: [1 x 1] integer, index of reference channel, e.g. 257
% MontageInfo: structure with channels coordinates, labels, etc. see for
%              example get_EGI_257ch_sensors_info.m
%
%  Outputs
% ---------
% Bad_Channels: [n x 1] indices of bad channels
% Good_Channels: [n x 1] indices of good channels
% Bad_Segments: [n x 1] indices of bad channels
% Good_Segments: [n x 1] indices of good channels
% Metrics: structure array containing metrics used for detecting bad channels &
%          segments, with corresponding values, and other information. Indices
%          of channels which were preliminarily rejected (e.g. flat channels) are
%          given a value of NaN for "Correlations", "Variance" and "Hurst_exp" metrics.
% Results: structure with some statistics on final detection results
% Settings: structure with "z"-score thresholds, window length & overlap
%           for bad segments detection.
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, May 2019
% NB: - this function focuses on steps 1 & 2 of FASTER pipeline (i.e. no
%       ICA is performed within this function
%     - input EEG traces are assumed to be already filtered
%     - modification compared to official FASTER pipeline: here, we do not
%       make epochs yet, so we apply step 2 of FASTER (bad epochs rejection,
%       which is necessary before applying ICA) using a sliding-window approach:
%       EEG traces are epoched in highly overlapping segments of 2 seconds, and
%       we reject segments for which metrics are bad in all overlapping time
%       points
%     - Compared to the original method, I found that the measures and
%       threshold used in the original study are not appropriate for my
%       dataset. This could be due to the fact that the data come from
%       epileptic patients... ? I changed the weights for each metric,
%       because the variance looks accurate to find bad channels, while
%       other measures, especially the Hurst exponent, give a lot of false
%       positives.
%     - In addition, compared to the original method, this function is
%       called recursively: this is because bad epochs will also influence
%       bad channels detection by biasing the metrics.
%
%       The principle is thus to do something as follows:
%
%                   ,- > bad channels detection & removal _
%                  /                                        \
%                 |                                          |
%                 |                                         /
%                  \_  bad segments detection & removal  <-´
%
%       My pipeline for bad channels and bad segments detection therefore
%       looks as follows:
%
%              STAGE 1
%           => bad channels detection & removal
%               => bad segments detection and removal (without bad channels)
%                      STAGE 2
%                   => bad channels re-detection without bad segments
%                       => bad segments re-detection without updated bad channels
%
%       Moreover: a) this implementation adds preliminary bad segments
%       detection using GFP (Global Field Power) and Dissimilarity, as well
%       as preliminary bad channels (without segments preliminarily
%       detected as bad) using the variance across channels in epochs: if
%       any channel is above a certain z-score, the whole epoch is
%       rejected; b) a new algorithm for detection of segments with single (or
%       few) channel(s) outlier signal is used. The principle of the latter
%       is desperately simple: it calculates, at every time point, the
%       maximal difference (diff) between SORTED potentials. The metric is
%       normalized using robust z-score, and if it is above a critical
%       value, the epoch is rejected.
%     - Uses now robust z-scores ((X - median(X))/mad(X)) instead of z-scores ((X - mean(X))/std(X))
%     - Uses average absolute difference across time for each channel to
%       detect noisy channels. Leave-one-out GFP and dissimilarity were
%       promising metrics but detect also good electrodes, because it looks
%       at channels which contribute to average GFP and dissimilarity
%       values across the whole recording, so this method was dropped.
%     - In some cases, there are extremely extreme values, e.g. around 200'000,
%       leading to robust z-scores above 100. In this case, huge artefacts prevent
%       the detection of small artefacts... I added a recursive part for
%       this section, and the indices of bad time frames update and feeds
%       itself such that until the maximal z-score is above 100, it
%       recalculate MaxDiff on remaining time frames. This technique could
%       potentially be used for other metrics, but so far I did not observe
%       measures with a distribution as skewed as this one in some cases,
%       so it might also just add unnecessary overhead and computational
%       complexity.
%-------------------------------------------------------------------------
% #TODO, #RM@FBMlab: include entropy measure using Murphy's code but scale
%                    before: entropy(scale_0_1(EEG(AllChansButRef,:)'),1)
% #TODO, #RM@FBMlab: if two bad segments are separated by less than 500 ms,
% consider the segment in between as bad as well... ?
%-------------------------------------------------------------------------

%=========================================
%========== Critical "z"-scores ==========
%=========================================
%% Channels:
Settings.Zcrit.Chan.Var = 15;%10;%3; % original threshold, because when there are bad channels it works well
% #RM@FBMlab, 2019-05-17: original z-threshold of 3 was changed back to a
% higher value, because "ok" channels were marked as bad! => from 3 to 5 =>
% from 5 to 10, because it marks good channels as bad! => with robust
% z-score, this metric sometimes pick channels that are good, raised the
% threshold up to 15!
Settings.Zcrit.Chan.Corr = 15;%10; % this measure does not work very well anyhow => this metric usually does not really help, raised threshold to 15!
Settings.Zcrit.Chan.Hurst = 15;%10; % this measure does not work very well anyhow => with robust
% z-score, this metric sometimes pick channels that are good, raised the
% threshold up to 15!
Settings.Zcrit.Chan.Std = 10;%5;%10; % changed it to 3 because bad channels were not picked with too conservative threshold.
% Anyhow, like for GFP and Diss, this might not really refer to a z-score
% anymore, because we restore baseline to 0 after having applied robust
% z-score and z-score...
% Settings.Zcrit.Chan.LeaveOneOutGFP = 10;
% Settings.Zcrit.Chan.LeaveOneOutDiss = 10;
Settings.Zcrit.Chan.AvgAbsDiff = 4;

%% Epochs:
Settings.Zcrit.Ep.Dev = 10;%5;
Settings.Zcrit.Ep.Var = 5;%10;%5; => missing some really bad segments otherwise, that neither MaxDiff nor GFP are catching !
Settings.Zcrit.Ep.AmpDiff = 10;%5;%3;%5; % 2019-05-16: testing lower threshold to see if false negatives disappear in segments without abrupt changes but large voltage differences in many electrodes, did not change the results so much...
Settings.Zcrit.Ep.GFP = 10;%8;%10;
Settings.Zcrit.Ep.Diss = 10;%8;%10;%9;%10;

%% Channels in epochs:
% If one has many channels and only few channels are bad, the weight of
% all other channels might prevent detection of bad segments... Use the
% same principle but apply them for each electrode at every epoch, if any
% electrode shows as bad, exclude it.
Settings.Zcrit.EpChan.Var = 5;%10; => missing some really bad segments otherwise, that neither MaxDiff nor GFP are catching !
% Zcrit.EpChan.MaxDevMean = 6;

% My method (maximal difference between sorted voltages across electrodes
% for each time point):
Settings.Zcrit.EpChan.MaxDiff = 10;%7;%15;%3;;%5; % => with robust z-score, this metric was marking most time frames as bad, because the threshold was set to a really low value (1.5 or 2), setting it to 5 is enough now that we switched from standard z-score to robust z-score
Settings.Zcrit.EpChan.MaxDiffMaxIter = 100;
Settings.Zcrit.EpChan.MaxDiffMax = 20;

% Moving-segments detection
Settings.SegmentWindowLengthRange = [0.8 1]; % in seconds, initially 1.5 to 2
Settings.SegmentsOverlap = 0.95; % desired level of overlap between segments

%=====================================================
%% Iterative bad channels / segments detection process
%=====================================================

% For what we aim to achieve here, single precision is largely
% sufficient, whereas double precision will just increase computational
% load and time:
EEG = single(EEG);

tic;
fprintf('\n=============== Iteration 1/2... ===============\n')

%% Preliminary bad segments rejection
% This avoids to bias bad channels detection, because once channels are
% out, the whole process is biased and it is more difficult to go back.
% Here we want something conservative to remove only really bad segments.
% This is refactored FBMlab code, to which I added the Dissimilarity metric,
% because it can pick bad epochs that are not detected by GFP.
% Segments spotted by this step are usually extremely bad, will
% drastically change the results of the bad channels detection and should
% be removed.

fprintf('Preliminary bad segments detection (GFP & dissimilarity)...\n')
[BadSeg,GFP,Diss] = GFP_Diss_bad_seg(EEG,SamplingFreq,Settings.Zcrit.Ep.GFP,Settings.Zcrit.Ep.Diss);

fprintf('Found %.2f %% of TF...\n',sum(BadSeg)/length(BadSeg)*100)

%% Segments with outlier electrodes
% This is a metric I suggest to detect segments in which one electrode
% among all goes crazy, while others stay normal. These segments would
% typically be considered bad when doing manual marking, but neither the
% GFP, nor the Dissimilarity will catch these, because the more electrodes
% one has, the less one electrode will weight in the result of the GFP and
% Dissimilarity. With high-density EEG, the proposed metric (maximal
% difference among electrodes voltage) will catch segments with outlier
% electrode behaviour particularly well.

fprintf('Maximal amplitude difference between sorted potentials (recursive)...\n')
MaxDiff = nan(size(EEG,2),1);
for t = 1:size(EEG,2)
    MaxDiff(t) = max(diff(sort(EEG(:,t))));
end
ScaledMaxDiff = ones(size(EEG,2),1)*1e3;
Count = 1;
SegWithOutChanTemp = false(size(EEG,2),1);
while max(ScaledMaxDiff(:))>Settings.Zcrit.EpChan.MaxDiffMax && Count<Settings.Zcrit.EpChan.MaxDiffMaxIter
    fprintf('Iteration %d (max = %d)...\n',Count,Settings.Zcrit.EpChan.MaxDiffMaxIter)
    EEGtemp = EEG(:,~SegWithOutChanTemp);
%     size(EEGtemp,2)
    MaxDiffTemp = nan(sum(~SegWithOutChanTemp),1);
    for t = 1:size(EEGtemp,2)
        MaxDiffTemp(t) = max(diff(sort(EEGtemp(:,t))));
    end
    Count = Count+1;
    ScaledMaxDiff = MAD_normalize_single_by_col(MaxDiffTemp);
    SegWithOutChanTemp(~SegWithOutChanTemp) = (SegWithOutChanTemp(~SegWithOutChanTemp) + ScaledMaxDiff>Settings.Zcrit.EpChan.MaxDiff)>0;
%     sum(SegWithOutChanTemp)
end
SegWithOutChan = smooth(single(SegWithOutChanTemp(:)),round(SamplingFreq/10))>0;
BadSeg(SegWithOutChan)=true;

fprintf('Found %.2f %% of TF...\n',sum(SegWithOutChan)/length(SegWithOutChan)*100)

%% Checkpoint: check bad segments:
% If almost all segments were rejected, it means that there are few (or
% even just one) electrodes driving the results of bad segments detection,
% so we forget all bad segments and goes on with bad channels detection. We
% will catch bad segments again later...
if sum(BadSeg)>(0.95*length(BadSeg))
    warning('Wow, >95% of EEG time points were detected as bad... You might have some real crazy channel throughout the whole recording, ignoring preliminary bad segments for now...')
    BadSeg = false(size(BadSeg));
end

%% Preliminary bad channels detection
% Again, this step avoids to bias subsequent steps. Here we mostly want to
% exclude flat and really jerky channels that would affect next metrics.
% Based on FBMlab code, adapted threshold to remove only extremely bad
% channels.
% The channels detected here are either flat, broken or extremely noisy.
% Compared to original FBM lab code, we do not include the reference
% channel in the calculation, because it would always be flagged as bad (all 0).

fprintf('Preliminary bad channels detection (std)...\n')

% Flat and jerky channels:
[FlatJerkyE,std_ch] = std_bad_chans(EEG,RefChan,Settings.Zcrit.Chan.Std);

fprintf('Found %d bad channels...\n',sum(FlatJerkyE))

%% Bad channels detection
% Here we estimate the average GFP and dissimilarity in a
% leave-one-out fashion, such that channels with outlier behaviour are
% spotted. These channels might look like EEG but have very low informational quality.
% They stand out from the others, even after filtering, and therefore
% greatly degrade the GFP (and dissimilarity).
% In addition, we use some of FASTER metrics to spot bad channels, but are
% rather conservative in our method because we want to avoid excluding good
% channels: in our experience, FASTER metrics often consider good channels
% as bad, so we use a high threshold and reject channels only if they
% exhibit extremely unlikely FASTER metrics values.

fprintf('Detecting bad channels (mean absolute derivative)...\n')

% Leave-one-out Dissimilarity and GFP => this method sometimes picks good channels
% => some channels, when removed, will indeed change the GFP but not in a good way!
% [NoisyChans,L1Ogfp,L1Odiss] = L1O_GFP_Diss(EEG,RefChan,BadSeg,...
%     Settings.Zcrit.Chan.LeaveOneOutGFP,Settings.Zcrit.Chan.LeaveOneOutDiss);
AvgAbsDiff = mean(abs(diff(EEG')));
NoisyChans = MAD_normalize_single_by_col(vect(AvgAbsDiff))>Settings.Zcrit.Chan.AvgAbsDiff;

fprintf('Found %d channels...\n',sum(NoisyChans))

% FASTER metrics (only on channels that passed the previous tests,
% otherwise this will bias the measures):
[CorrelationsTemp, VarianceTemp, Hurst_expTemp,eeg_chans,...
    s_pol_dist,dist_inds,idist_inds, proj_dist, xyz_dist, pol_dist] = ...
    metrics_chans(EEG(:,~BadSeg),RefChan, (FlatJerkyE+NoisyChans)>0, MontageInfo);

MetricsChan = abs(MAD_normalize_single_by_col([CorrelationsTemp(:),VarianceTemp(:),Hurst_expTemp(:)]));
% to look also at Z < threshold, not only Z > threshold

ChanThresh = [repmat(Settings.Zcrit.Chan.Corr,size(MetricsChan,1),1),...
    repmat(Settings.Zcrit.Chan.Var,size(MetricsChan,1),1),...
    repmat(Settings.Zcrit.Chan.Hurst,size(MetricsChan,1),1)];

Bad_Chans = (FlatJerkyE+NoisyChans)>0;
% Among channels that are not really bad, which are identified as bad by
% FASTER metrics:
BadFASTER = sum(MetricsChan>ChanThresh,2)>0;
Bad_Chans(~Bad_Chans) = BadFASTER;

fprintf('Found %d channels...\n',sum(BadFASTER))

if ~isempty(intersect(find(Bad_Chans),RefChan)) && ~isempty(RefChan)
    error('Reference channel was detected as bad but it should not be the case, please check for issues in the code!')
end

%% Checkpoint: check bad channels:
% If almost all channels were rejected, that likely means we missed some
% bad segments, so we forget bad channels and will catch them again
% later...
if sum(Bad_Chans)>(size(EEG,1)-2)
    warning('Wow, (almost) all channels look bad... You might have a lot of bad epochs throughout the whole recording, ignoring preliminary bad channels for now...')
    Bad_Chans = false(size(Bad_Chans));
end

Bad_Channels = find(Bad_Chans);
Good_Channels = setdiff(eeg_chans,Bad_Channels)';

Correlations = nan(size(EEG,1),1);
Variance = nan(size(EEG,1),1);
Hurst_exp = nan(size(EEG,1),1);
% This way we avoid wrong indices for metrics...
Correlations((FlatJerkyE+NoisyChans)<1) = CorrelationsTemp;
Variance((FlatJerkyE+NoisyChans)<1) = VarianceTemp;
Hurst_exp((FlatJerkyE+NoisyChans)<1) = Hurst_expTemp;

Metrics.Correlations = Correlations;
Metrics.Variance = Variance;
Metrics.Hurst_exp = Hurst_exp;
Metrics.pol_dist = pol_dist;
Metrics.xyz_dist = xyz_dist;
Metrics.proj_dist = proj_dist;
Metrics.GFP = GFP;
Metrics.Dissimilarity = Diss;
Metrics.std_ch = std_ch;
Metrics.MaxDiff = MaxDiff;
% Metrics.L1Ogfp = L1Ogfp;
% Metrics.L1Odiss = L1Odiss;
Metrics.AvgAbsDiff = AvgAbsDiff;

%% Bad segments detection
% Here we use some of FASTER metrics to detect bad epochs, but we apply it
% in a moving-window fashion (with substantial overlap between segments).
% This way, we do not need epoching before next processing stages and we
% can mark all continuous EEG. If epochs change later, we won't have to
% re-do this part. However, we use very conservative threshold for FASTER
% metrics, because, again, in our experience, the default values (even with
% standard z-scores instead of robust z-scores) will reject too many good
% segments (such as the ones with strong alpha waves: if we want to reject
% these at all, it is not the aim of this algorithm, but rather what ICA
% would do).
% In addition, we use FASTER's variance metric (which is originally used
% for bad channels detection) to each segment, and we consider a segment to
% be bad if any channel shows abnormal variance. We use the same
% moving-window approach as above. This is heavy but works relatively well
% (especially compared to other metrics) and is reasonably fast to compute.

fprintf('Detecting bad segments (deviation, variance, amplitude range)...\n')

[Deviations, EpVariance, AmpDiffs, Dev2samp, EpVar2samp, AmpDiff2samp, CritThresh] = ...
    metrics_epochs(EEG(Good_Channels,:),Settings.SegmentWindowLengthRange,SamplingFreq,Settings.SegmentsOverlap,Settings.Zcrit);

BadFASTER = any([Dev2samp,EpVar2samp,AmpDiff2samp]>=repmat(CritThresh,1,3),2);
Bad_Segments = (BadSeg+BadFASTER)>0;

fprintf('Found %.2f %% of TF...\n',sum(BadFASTER)/length(BadFASTER)*100)

fprintf('Detecting segments with bad channels (variance)...\n')

[ChanInEpVariance, Var2samp, CritThresh] = ...
    metrics_epochs_chans(EEG(Good_Channels,:),Settings.SegmentWindowLengthRange,SamplingFreq,...
    Settings.SegmentsOverlap,Settings.Zcrit,RefChan,s_pol_dist,dist_inds,idist_inds);

BadFASTER = any(Var2samp>=CritThresh,2);

Bad_Segments = (Bad_Segments+BadFASTER)>0;
% Good_Segments = ~Bad_Segments;

fprintf('Found %.2f %% of TF...\n',sum(BadFASTER)/length(BadFASTER)*100)

Metrics.Deviations = Deviations;
Metrics.EpVariance = EpVariance;
Metrics.AmpDiffs = AmpDiffs;
Metrics.Dev2samp = Dev2samp;
Metrics.EpVar2samp = EpVar2samp;
Metrics.AmpDiff2samp = AmpDiff2samp;
Metrics.CritThresh = CritThresh;
Metrics.ChanInEpVariance = ChanInEpVariance;
Metrics.Var2samp = Var2samp;

%% SECOND STAGE: update bad channels and bad segments
fprintf('\n=============== Iteration 2/2... ===============\n')
%% Second-stage bad channels detection
% Now that bad segments were excluded, we re-detect bad channels with
% updated EEG traces: channels that were excluded earlier, just because of
% multiple bad segments, might come back.
% Bad segments detection using GFP, dissimilarity, MaxDiff & std is moved to after bad
% channels detection in the 2nd stage:
fprintf('Re-detecting bad channels (std)...\n')
[FlatJerkyE,std_ch] = std_bad_chans(EEG(:,~Bad_Segments),RefChan,Settings.Zcrit.Chan.Std);
fprintf('Found %d bad channels...\n',sum(FlatJerkyE))
fprintf('Re-detecting bad channels (mean absolute derivative)...\n')
AvgAbsDiff = mean(abs(diff(EEG(:,~Bad_Segments)')));
NoisyChans = MAD_normalize_single_by_col(vect(AvgAbsDiff))>Settings.Zcrit.Chan.AvgAbsDiff;
fprintf('Found %d channels...\n',sum(NoisyChans))
[CorrelationsTemp, VarianceTemp, Hurst_expTemp,eeg_chans,...
    s_pol_dist,dist_inds,idist_inds, proj_dist, xyz_dist, pol_dist] = ...
    metrics_chans(EEG(:,~Bad_Segments),RefChan, (FlatJerkyE+NoisyChans)>0, MontageInfo);
MetricsChan = abs(MAD_normalize_single_by_col([CorrelationsTemp(:),VarianceTemp(:),Hurst_expTemp(:)]));
ChanThresh = [repmat(Settings.Zcrit.Chan.Corr,size(MetricsChan,1),1),...
    repmat(Settings.Zcrit.Chan.Var,size(MetricsChan,1),1),...
    repmat(Settings.Zcrit.Chan.Hurst,size(MetricsChan,1),1)];
Bad_Chans = (FlatJerkyE+NoisyChans)>0;
BadFASTER = sum(MetricsChan>ChanThresh,2)>0;
Bad_Chans(~Bad_Chans) = BadFASTER;
fprintf('Found %d channels...\n',sum(BadFASTER))
if ~isempty(intersect(find(Bad_Chans),RefChan)) && ~isempty(RefChan)
    error('Reference channel was (re-)detected as bad but it should not be the case, please check for issues in the code!')
end
if sum(Bad_Chans)>(size(EEG,1)-2)
    warning('Wow, (almost) all channels look bad... You might have a lot of bad epochs throughout the whole recording, ignoring preliminary bad channels for now...')
    Bad_Chans = false(size(Bad_Chans));
end
Bad_Channels = find(Bad_Chans);
Good_Channels = setdiff(eeg_chans,Bad_Channels)';
Correlations = nan(size(EEG,1),1);
Variance = nan(size(EEG,1),1);
Hurst_exp = nan(size(EEG,1),1);
Correlations((FlatJerkyE+NoisyChans)<1) = CorrelationsTemp;
Variance((FlatJerkyE+NoisyChans)<1) = VarianceTemp;
Hurst_exp((FlatJerkyE+NoisyChans)<1) = Hurst_expTemp;
Metrics.Correlations = Correlations;
Metrics.Variance = Variance;
Metrics.Hurst_exp = Hurst_exp;
Metrics.pol_dist = pol_dist;
Metrics.xyz_dist = xyz_dist;
Metrics.proj_dist = proj_dist;
Metrics.std_ch = std_ch;
% Metrics.L1Ogfp = L1Ogfp;
% Metrics.L1Odiss = L1Odiss;
Metrics.AvgAbsDiff = AvgAbsDiff;

%% Second-stage bad segments detection
% Now that we have updated bad channels, we re-detect bad segments, given
% that remaining channels should be good apart from some segments
% previously detected. We should find again more or less the same segments,
% but possibly slightly more.
fprintf('Re-detecting bad segments (GFP & dissimilarity)...\n')
[BadSeg,GFP,Diss] = GFP_Diss_bad_seg(EEG(Good_Channels,:),SamplingFreq,Settings.Zcrit.Ep.GFP,Settings.Zcrit.Ep.Diss);
fprintf('Found %.2f %% of TF...\n',sum(BadSeg)/length(BadSeg)*100)

fprintf('Maximal amplitude difference between sorted potentials (recursive)...\n')
MaxDiff = nan(size(EEG(Good_Channels,:),2),1);
for t = 1:size(EEG(Good_Channels,:),2)
    MaxDiff(t) = max(diff(sort(EEG(Good_Channels,t))));
end
ScaledMaxDiff = ones(size(EEG,2),1)*1e3;
Count = 1;
SegWithOutChanTemp = false(size(EEG,2),1);
while max(ScaledMaxDiff(:))>Settings.Zcrit.EpChan.MaxDiffMax && Count<Settings.Zcrit.EpChan.MaxDiffMaxIter
    fprintf('Iteration %d (max = %d)...\n',Count,Settings.Zcrit.EpChan.MaxDiffMaxIter)
    EEGtemp = EEG(Good_Channels,~SegWithOutChanTemp);
%     size(EEGtemp,2)
    MaxDiffTemp = nan(sum(~SegWithOutChanTemp),1);
    for t = 1:size(EEGtemp,2)
        MaxDiffTemp(t) = max(diff(sort(EEGtemp(:,t))));
    end
    Count = Count+1;
    ScaledMaxDiff = MAD_normalize_single_by_col(MaxDiffTemp);
    SegWithOutChanTemp(~SegWithOutChanTemp) = (SegWithOutChanTemp(~SegWithOutChanTemp) + ScaledMaxDiff>Settings.Zcrit.EpChan.MaxDiff)>0;
%     sum(SegWithOutChanTemp)
end
SegWithOutChan = smooth(single(SegWithOutChanTemp(:)),round(SamplingFreq/10))>0;
fprintf('Found %.2f %% of TF...\n',sum(SegWithOutChan)/length(SegWithOutChan)*100)
BadSeg(SegWithOutChan)=true;
fprintf('Detecting bad segments (deviation, variance, amplitude range)...\n')
[Deviations, EpVariance, AmpDiffs, Dev2samp, EpVar2samp, AmpDiff2samp, CritThresh] = ...
    metrics_epochs(EEG(Good_Channels,:),Settings.SegmentWindowLengthRange,SamplingFreq,Settings.SegmentsOverlap,Settings.Zcrit);
BadFASTER = any([Dev2samp,EpVar2samp,AmpDiff2samp]>=repmat(CritThresh,1,3),2);
Bad_Segments = (BadSeg+BadFASTER)>0;
fprintf('Found %.2f %% of TF...\n',sum(BadFASTER)/length(BadFASTER)*100)
fprintf('Re-detecting segments with bad channels (variance)...\n')
[ChanInEpVariance, Var2samp, CritThresh] = ...
    metrics_epochs_chans(EEG(Good_Channels,:),Settings.SegmentWindowLengthRange,SamplingFreq,...
    Settings.SegmentsOverlap,Settings.Zcrit,RefChan,s_pol_dist,dist_inds,idist_inds);
BadFASTER = any(Var2samp>=CritThresh,2);
fprintf('Found %.2f %% of TF...\n',sum(BadFASTER)/length(BadFASTER)*100)
Bad_Segments = (Bad_Segments+BadFASTER)>0;
Good_Segments = ~Bad_Segments;
Metrics.Deviations = Deviations;
Metrics.EpVariance = EpVariance;
Metrics.AmpDiffs = AmpDiffs;
Metrics.Dev2samp = Dev2samp;
Metrics.EpVar2samp = EpVar2samp;
Metrics.AmpDiff2samp = AmpDiff2samp;
Metrics.CritThresh = CritThresh;
Metrics.ChanInEpVariance = ChanInEpVariance;
Metrics.Var2samp = Var2samp;
Metrics.GFP = GFP;
Metrics.Dissimilarity = Diss;
Metrics.MaxDiff = MaxDiff;

% NB: we could continue this until we reach a steady state, but each
% detection is relatively time-consuming, and it is difficult to find a
% good trade-off between learning rate and stability. Thus, at this point,
% we stop and assume that most non-physiological artefacts are detected
% (and can be excised before performing ICA, which will remove the
% remaining artefacts (also physiological ones)).

%% Stats and other useful variables

% Construct onsets and offsets of bad segments:
OnsetsBadIdx = find(diff(Bad_Segments)>0);
OffsetsBadIdx = find(diff(Bad_Segments)<0);
% Correct for BOF & EOF:
if length(OnsetsBadIdx)>length(OffsetsBadIdx)
    OffsetsBadIdx(end+1) = size(EEG,2);
end
if length(OffsetsBadIdx)>length(OnsetsBadIdx)
    OnsetsBadIdxBKP = OnsetsBadIdx;
    OnsetsBadIdx = nan(size(OnsetsBadIdx)+[1 0]);
    OnsetsBadIdx(2:end) = OnsetsBadIdxBKP;
    OnsetsBadIdx(1) = 1;
end

fprintf('======================\nDetected bad segments:\n\n')
fprintf('\tFrom:\t\tTo:\n')
disp([OnsetsBadIdx,OffsetsBadIdx])

fprintf('======================\n%.2f %% of time points are bad.\n\n',sum(Bad_Segments)/length(Bad_Segments)*100)

fprintf('======================\n%d bad channels detected:\n\n',length(Bad_Channels))
disp(MontageInfo.clinicalname(Bad_Channels))

Results.Segments.Onsets = OnsetsBadIdx;
Results.Segments.Offsets = OffsetsBadIdx;
Results.Segments.NbadTF = sum(Bad_Segments);
Results.Segments.NgoodTF = sum(Good_Segments);
Results.Segments.PercentBadTF = sum(Bad_Segments)/length(Bad_Segments);

Results.Channels.Nbad = length(Bad_Channels);
Results.Channels.NamesBad = MontageInfo.clinicalname(Bad_Channels);
Results.Channels.NamesGood = MontageInfo.clinicalname(Good_Channels);
Results.Channels.PercentBad = length(Bad_Channels)/size(EEG,1);

Toc = toc;
if Toc<60
    fprintf('Bad channels and segments detection\nperformed in %d seconds.\n',round(Toc));
elseif (Toc/60)<60
    fprintf('Bad channels and segments detection\nperformed in %d minutes.\n',round(Toc/60));
else
    fprintf('Bad channels and segments detection\nperformed in %d hours and %d minutes.\n',round(Toc/3600),round(rem(Toc,3600)/60));
end
end

%% LOCAL-FUNCTIONS

%% SUB-FUNCTION FOR REALLY BAD SEGMENTS DETECTION
function [BadSeg,GFP,Diss] = GFP_Diss_bad_seg(EEG,SamplingFreq,ThreshGFP,ThreshDiss)
% 2 differences with the original method from FBMlab: uses robust z-score
% instead of z-score (median instead of mean, median absolute deviation
% instead of std); subtract the minimum from the scaled measure to restore
% back to 0, otherwise we get baseline below 0 because distribution is not
% normal at all!

[Diss,~,GFP] = get_Dissimilarity(EEG);
ScaledGFP = MAD_normalize_single_by_col(GFP(:));
ScaledDiss = MAD_normalize_single_by_col(Diss(:));

% Alternative, 2019-05-16: 
% The distribution of these measures is not normal at all, so
% restore it to 0 with -min:
% ScaledGFP = ScaledGFP-min(ScaledGFP);
% ScaledDiss = ScaledDiss-min(ScaledDiss);
% NB: the above is ok because there should not be negative values in GFP
% and Diss by definition.

Bad1 = ScaledGFP>ThreshGFP; % 10 is safe, otherwise we might run into false positives
Bad2 = ScaledDiss>ThreshDiss;
BadSeg = (Bad1+Bad2)>0;
% Let's be conservative and reject data close to the detected periods,
% because the indices usually are scattered around artifacts (this is
% because GFP and Diss will "oscillate" around them). We will thus remove
% also time points around it using a window of 1 second:
BadSeg = smooth(single(BadSeg(:)),SamplingFreq)>0;

end

%% SUB-FUCTION FOR FLAT & JERKY CHANNELS CHANNELS DETECTION
function [BadE,std_ch] = std_bad_chans(EEG,RefChan,Thresh)
% 2 differences with the original method from FBMlab: uses robust z-score
% instead of z-score (median instead of mean, median absolute deviation
% instead of std); uses a z-score threshold of 10 instead of 3 (the latter
% gives too much false positives, we just to exlude channels that are really crazy)

std_ch=max(abs(EEG(setdiff(1:size(EEG,1),RefChan),:)),[],2);
std_ch(~isfinite(std_ch))=0;
Bad1 = std_ch<eps('single'); % practically equals 0: we assume the data is at least float with single precision
ScaledStdCh = MAD_normalize_single_by_col(std_ch(:));
% ScaledStdCh = ScaledStdCh-min(ScaledStdCh);
Bad2 = ScaledStdCh>Thresh; % here the threshold is 10 instead of 3, we prefer to keep channels and remove epochs instead, and later, remove them if they are really bad.
BadE = (Bad1+Bad2)>0;
if ~isempty(RefChan)
    BadE(end+1)=false;
end

end

%% SUB-FUNCTION FOR NOISY AND OUTLIER CHANNELS DETECTION:

% % The following seemed to work well but catches sometimes good
% % channels...
% function [BadE,L1OGFP,L1Odiss] = L1O_GFP_Diss(EEG,RefChan,BadSeg,Thresh1,Thresh2)
% % Probably what Cartool does when it computes "a lot of stats"... This aims
% % to detect noisy / outlier channels...
% 
% fprintf('Leave-1-out GFP & dissimilarity...\n')
% eeg_chans = setdiff(1:size(EEG,1),RefChan);
% [L1Odiss,L1OGFP] = deal(nan(length(eeg_chans),1));
% for c = 1:length(eeg_chans)
%     prc_for_loop(c,length(eeg_chans),1);
%     L1o = setdiff(eeg_chans,c);
%     % Here in particular, single precision will be crucial: it allows a
%     % 25-40% decrease in computational time, whereas double precision will
%     % just increase computational load and time (each iteration takes 2-2.3
%     % seconds on my machine with double precision, and this will be done
%     % for each channel...)
%     [Diss,~,GFP] = get_Dissimilarity(EEG(L1o,~BadSeg));
%     L1OGFP(c) = mean(GFP);
%     L1Odiss(c) = mean(Diss);
% end
% 
% BadE = abs(MAD_normalize_single_by_col(L1OGFP(:)))>Thresh1;
% BadE = (BadE + (abs(MAD_normalize_single_by_col(L1Odiss(:))) > Thresh2) ) > 0;
% 
% if ~isempty(RefChan)
%     BadE(end+1)=false;
% end
% 
% end

%% SUB-FUNCTION FOR BAD-CHANNELS DETECTION:
function [Correlations, Variance, Hurst_exp,eeg_chans,s_pol_dist,dist_inds,idist_inds, proj_dist, xyz_dist, pol_dist] = metrics_chans(EEG,RefChan, BadE, EGI)

%% Weighting of metric by distance
eeg_chans = find(~BadE);

if ~isempty(RefChan) && length(RefChan)==1
    [pol_dist, xyz_dist, proj_dist]=elec_dist_mat_pol(EGI);
    [s_pol_dist, dist_inds] = sort(pol_dist(RefChan,eeg_chans));
    [~, idist_inds] = sort(dist_inds);
else
    [pol_dist,s_pol_dist,xyz_dist,proj_dist,dist_inds,idist_inds] = deal([]);
end

% Calculate correlations
% calc_indices=setdiff(eeg_chans,RefChan);
% ignore_indices=intersect(eeg_chans,RefChan);
if ~isempty(RefChan)
    ignore_indices = eeg_chans==RefChan;
else
    ignore_indices = [];
end
% corrs = abs(corrcoef(EEG(setdiff(eeg_chans,RefChan),:)'));
% mcorrs=zeros(size(eeg_chans));
% for u=1:length(calc_indices)
%     mcorrs(calc_indices(u))=mean(corrs(u,:));
% end
% mcorrs(ignore_indices)=mean(mcorrs(calc_indices));
corrs = abs(corrcoef(EEG(eeg_chans,:)'));
mcorrs = nanmean(corrs,2);
mcorrs(isnan(mcorrs))=nanmean(mcorrs);

% Quadratic correction for distance from reference electrode
if (~isempty(RefChan) && length(RefChan)==1)
    p = polyfit(vect(s_pol_dist),vect(mcorrs(dist_inds)),2);
    fitcurve = polyval(p,s_pol_dist);
    Correlations = vect(mcorrs(dist_inds)) - vect(fitcurve(idist_inds));
else
    Correlations = mcorrs;
end
Correlations(ignore_indices)=mean(Correlations); % #RM@FBMlab: otherwise the reference channel might be picked as bad because of the weighting...

% 3 Variance of the channels
vars = var(EEG(eeg_chans,:)'); %#ok<UDIM>
vars(~isfinite(vars))=mean(vars(isfinite(vars)));

% Quadratic correction for distance from reference electrode
if (~isempty(RefChan) && length(RefChan)==1)
    p = polyfit(vect(s_pol_dist),vect(vars(dist_inds)),2);
    fitcurve = polyval(p,s_pol_dist);
    corrected = vect(vars) - vect(fitcurve(idist_inds));
    
    Variance = corrected;
else
    Variance = vars;
end

Hurst_exp = nan(length(eeg_chans),1);
% 4 Hurst exponent
for u=1:length(eeg_chans)
    Hurst_exp(u) = hurst_exponent(EEG(eeg_chans(u),:));
end

% Correlations(isnan(Correlations)) = nanmean(Correlations);
% Variance(isnan(Variance)) = nanmean(Variance);
% Hurst_exp(isnan(Hurst_exp)) = nanmean(Hurst_exp);
% Correlations = Correlations-median(Correlations);
% Variance = Variance-median(Variance);
% Hurst_exp = Hurst_exp-median(Hurst_exp);

end

%% SUB-FUNCTION FOR BAD SEGMENTS DETECTION:
function [Deviations, EpVariance, AmpDiffs, Dev2samp, EpVar2samp, AmpDiff2samp, CritThresh] = metrics_epochs(EEG,SegmentWindowLengthRange,SamplingFreq,SegmentsOverlap,Zcrit)

Nsamp = size(EEG,2);
% Determine optimal segment duration:
SegmentsDurationsStepsWindow = ...
    unique(round((round(SamplingFreq*SegmentWindowLengthRange(1)):round(SamplingFreq*SegmentWindowLengthRange(2)))*(1-SegmentsOverlap))); % between 1.5 and 2 seconds
%     unique(round((round(SamplingFreq*1.5):round(SamplingFreq*2))*(1-SegmentsOverlap))); % between 1.5 and 2 seconds
[Remainder,IdxDur] = min(rem(Nsamp,SegmentsDurationsStepsWindow));
if Remainder~=0
    [Remainder,IdxDur] = max(rem(Nsamp,SegmentsDurationsStepsWindow));
end
StepsDuration = SegmentsDurationsStepsWindow(IdxDur);
SegmentsDuration = round(StepsDuration/(1-SegmentsOverlap));

Nseg = (Nsamp-Remainder)/StepsDuration;
if mod(Nseg,1)>0
    error('Could not calculate the number of segments based on the optimal window size, please check code...')
end
Nseg = Nseg-(SegmentsDuration/StepsDuration-1)+1;

Means = mean(EEG,2);
Deviations = nan(Nseg,1);
EpVariance = nan(Nseg,1);
AmpDiffs = nan(Nseg,1);
Count = 1;
for s = 1:Nseg
    prc_for_loop(s,Nseg,100);
    if s == Nseg
        EEGseg = EEG(:,Count:end);
    else
        EEGseg = EEG(:,Count:Count+SegmentsDuration-1);
    end
    Deviations(s) = mean(abs(mean(EEGseg,2)-Means));
    EpVariance(s) = mean(var(EEGseg,0,2));
    AmpDiffs(s) = mean(max(EEGseg,[],2)-min(EEGseg,[],2));
    Count = Count+StepsDuration;
end

Deviations = MAD_normalize_single_by_col(Deviations(:));
EpVariance = MAD_normalize_single_by_col(EpVariance(:));
AmpDiffs = MAD_normalize_single_by_col(AmpDiffs(:));
Dev2samp = zeros(size(EEG,2),1,'single');
EpVar2samp = zeros(size(EEG,2),1,'single');
AmpDiff2samp = zeros(size(EEG,2),1,'single');
CritThresh = zeros(size(EEG,2),1,'single');
Count = 1;
for s = 1:Nseg
    %     prc_for_loop(s,Nseg,30);
    if s == Nseg
        Dev2samp(Count:end) = Dev2samp(Count:end)+repmat(Deviations(s)>Zcrit.Ep.Dev,Nsamp-Count+1,1);
        EpVar2samp(Count:end) = EpVar2samp(Count:end)+repmat(EpVariance(s)>Zcrit.Ep.Var,Nsamp-Count+1,1);
        AmpDiff2samp(Count:end) = AmpDiff2samp(Count:end)+repmat(AmpDiffs(s)>Zcrit.Ep.AmpDiff,Nsamp-Count+1,1);
        CritThresh(Count:end) = CritThresh(Count:end)+ones(Nsamp-Count+1,1);
    else
        Dev2samp(Count:Count+SegmentsDuration-1) = Dev2samp(Count:Count+SegmentsDuration-1)+repmat(Deviations(s)>Zcrit.Ep.Dev,SegmentsDuration,1);
        EpVar2samp(Count:Count+SegmentsDuration-1) = EpVar2samp(Count:Count+SegmentsDuration-1)+repmat(EpVariance(s)>Zcrit.Ep.Var,SegmentsDuration,1);
        AmpDiff2samp(Count:Count+SegmentsDuration-1) = AmpDiff2samp(Count:Count+SegmentsDuration-1)+repmat(AmpDiffs(s)>Zcrit.Ep.AmpDiff,SegmentsDuration,1);
        CritThresh(Count:Count+SegmentsDuration-1) = CritThresh(Count:Count+SegmentsDuration-1)+ones(SegmentsDuration,1);
    end
    Count = Count+StepsDuration;
end

end

%% SUB-FUNCTION FOR EPOCHS WITH BAD CHANNELS DETECTION:
function [Variance, Var2samp, CritThresh] = metrics_epochs_chans(EEG,SegmentWindowLengthRange,SamplingFreq,SegmentsOverlap,Zcrit,RefChan,s_pol_dist,dist_inds,idist_inds)

Nsamp = size(EEG,2);
% Determine optimal segment duration:
SegmentsDurationsStepsWindow = ...
    unique(round((round(SamplingFreq*SegmentWindowLengthRange(1)):round(SamplingFreq*SegmentWindowLengthRange(2)))*(1-SegmentsOverlap))); % between 1.5 and 2 seconds
%     unique(round((round(SamplingFreq*1.5):round(SamplingFreq*2))*(1-SegmentsOverlap))); % between 1.5 and 2 seconds
[Remainder,IdxDur] = min(rem(Nsamp,SegmentsDurationsStepsWindow));
if Remainder~=0
    [Remainder,IdxDur] = max(rem(Nsamp,SegmentsDurationsStepsWindow));
end
StepsDuration = SegmentsDurationsStepsWindow(IdxDur);
SegmentsDuration = round(StepsDuration/(1-SegmentsOverlap));

Nseg = (Nsamp-Remainder)/StepsDuration;
if mod(Nseg,1)>0
    error('Could not calculate the number of segments based on the optimal window size, please check code...')
end
Nseg = Nseg-(SegmentsDuration/StepsDuration-1)+1;

Variance = nan(Nseg,1);
% MaxDevMean = nan(Nseg,1);
Count = 1;
for s = 1:Nseg
    prc_for_loop(s,Nseg,50);
    if s == Nseg
        EEGseg = EEG(:,Count:end);
    else
        EEGseg = EEG(:,Count:Count+SegmentsDuration-1);
    end
    
    % Here we evaluate only variance, it is fast and works well:
    vars = var(EEGseg,[],2);
    vars(~isfinite(vars))=mean(vars(isfinite(vars)));
    
    % if bad electrodes were detected, indices for quadratic correction
    % based on distance from reference electrode must be adapted:
    ChanOK = length(vars);
    
    % Quadratic correction for distance from reference electrode
    if (~isempty(RefChan) && length(RefChan)==1)
        p = polyfit(vect(s_pol_dist(dist_inds<=ChanOK)),vect(vars(dist_inds(dist_inds<=ChanOK))),2);
        fitcurve = polyval(p,s_pol_dist(dist_inds<=ChanOK));
        corrected = vect(vars) - vect(fitcurve(idist_inds(idist_inds<=ChanOK)));
        
        VarTemp = corrected;
    else
        VarTemp = vars;
    end
    
    Variance(s) = max(VarTemp);
%     MeanEEG = mean(EEGseg,1);
%     MaxDevMeanTemp = nan(size(EEGseg,1),1);
%     MaxDevMeanTemp = corr(EEGseg',MeanEEG');
%     for c = 1:size(EEGseg,1)
%         MaxDevMeanTemp(c) = fastcorrelation(EEGseg(c,:),MeanEEG);
%     end
%     MaxDevMean(s) = nanmean(MaxDevMeanTemp);
    
    Count = Count+StepsDuration;
end

VarianceZ = MAD_normalize_single_by_col(Variance(:));
% MaxDevMeanZ = normalize_by_col(MaxDevMean-nanmedian(MaxDevMean));

Var2samp = zeros(size(EEG,2),1,'single');
% MaxDevMean2samp = zeros(size(EEG,2),1,'single');
CritThresh = zeros(size(EEG,2),1,'single');
Count = 1;
for s = 1:Nseg
    %     prc_for_loop(s,Nseg,30);
    if s == Nseg
        Var2samp(Count:end) = Var2samp(Count:end)+repmat(VarianceZ(s)>Zcrit.EpChan.Var,Nsamp-Count+1,1);
%         MaxDevMean2samp(Count:end) = MaxDevMean2samp(Count:end)+repmat(MaxDevMeanZ(s)>Zcrit.EpChan.MaxDevMean,Nsamp-Count+1,1);
        CritThresh(Count:end) = CritThresh(Count:end)+ones(Nsamp-Count+1,1);
    else
        Var2samp(Count:Count+SegmentsDuration-1) = Var2samp(Count:Count+SegmentsDuration-1)+repmat(VarianceZ(s)>Zcrit.EpChan.Var,SegmentsDuration,1);
%         MaxDevMean2samp(Count:Count+SegmentsDuration-1) = MaxDevMean2samp(Count:Count+SegmentsDuration-1)+repmat(MaxDevMeanZ(s)>Zcrit.EpChan.MaxDevMean,SegmentsDuration,1);
        CritThresh(Count:Count+SegmentsDuration-1) = CritThresh(Count:Count+SegmentsDuration-1)+ones(SegmentsDuration,1);
    end
    Count = Count+StepsDuration;
end

end



% [Bad_Chans, Good_Channels,...
%     ~, ~,...
%     Metrics(2)] = FASTER_improved(EEG(:,Good_Segments), SamplingFreq, RefChan, MontageInfo, Settings);
% 
% if ~isempty(RefChan)
%     RefChanIdx = find(Good_Channels==RefChan);
% else
%     RefChanIdx = [];
% end
% 
% [~, ~,...
%     Bad_Segments, Good_Segments,...
%     Metrics(3)] = FASTER_improved(EEG(Good_Channels,:), SamplingFreq, RefChanIdx, MontageInfo, Settings);

% fprintf('Iteration 4/5...\n')
% [Bad_Channels, Good_Channels,...
%     ~, ~,...
%     ~, Metrics(4)] = stripped_down_FASTER(EEG(:,Good_Segments), SamplingFreq, RefChan);
% 
% fprintf('Iteration 5/5...\n')
% [~, Good_Channels,...
%     Bad_Segments, Good_Segments,...
%     EGI, Metrics(5)] = stripped_down_FASTER(EEG(Good_Channels,:), SamplingFreq, RefChan);
