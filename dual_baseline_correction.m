function baseline_corrected_EEG = dual_baseline_correction( EEG, Trials, PreStim, PostStim, NotBaseline )
% dual_baseline_correction: baseline correction for EEG "subtracting EACH track with ITS average within a given time period" (Cartool reference guide)
%
% baseline_corrected_EEG = dual_baseline_correction( EEG, Trials, PreStim, PostStim, NotBaseline )
%
%  Inputs
% --------
% EEG         : channels x time frames EEG traces
% Trials      : trial timings (in samples)
% PreStim     : pre-stimulus duration (in samples)
% PostStim    : post-stimulus duration (in samples)
% NotBaseline : pre-stimulus duration to exclude for baseline-correction
%
%  Outputs
% ---------
% baseline_corrected_EEG
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, July 2018
%-------------------------------------------------------------------------

NumChan = size(EEG,1);
%======== baseline correction ========
% "subtracting EACH track with ITS average within a given time period" (Cartool reference guide)
baseline_corrected_EEG = nan(size(EEG,1),length(Trials),PreStim+PostStim+1); % number of electrodes x 20 epochs (stimulation pulses) x 1000 datapoints (ms)
for n = 1:NumChan
    for epoch = 1:length(Trials)
        baseline_corrected_EEG(n,epoch,:) = ...
            EEG(n,(Trials(epoch)-PreStim):(Trials(epoch)+PostStim))...
            -mean(EEG(n,(Trials(epoch)-PreStim):Trials(epoch)-NotBaseline));
    end
end

end

