function [SigEpochs, SigEpochsAmp] = sig_mat(EEG, HeightThresh, WidthTresh, PreStim, PostStim, NotBaseline)
% SIG_MAT : get significant epochs from EEG trace
% 
% Usage:
%-------------------------------------------------------------------------
% SigEpochs = sig_mat(EEG, HeightThresh, WidthTresh, PreStim, PostStim, NotBaseline)
%
%
% Inputs:
%-------------------------------------------------------------------------
% EEG: channels x time points array of EEG traces
%
% HeightTresh : number of std defining the baseline (one value that will be
%               multiplied by the std of each electrode)
%
% WidthTresh : minimal number of time points during which a period should be
%              above HeightThresh
%
% PreStim : number of time points before stimulation onset
%
% PostStim : number of time points after stimulation onset (should in
%            principle equal number of columns in EEG minus PreStim -1
%
% NotBaseline : number of samples during which the pre-stimulus period
%               cannot be considered as a baseline anymore (see Megevand et
%               al (https://doi.org/10.1089/brain.2017.0527)
%
%
% Outputs:
%-------------------------------------------------------------------------
% SigEpochs : PostStim + 1 x number of channels binary array of
%             significance
% SigEpochsAmp : same as SigEpochs except that if SigEpochs == 1 the value
% in SigEpochsAmp is equal to the difference between EEG and HeightTresh *
% std(EEG) during baseline
%
% Nota bene:
%-------------------------------------------------------------------------
% Baseline is defined per channel, like in:
% Keller et al, 2011 (http://doi.org/10.1073/pnas.1019750108)
%
% Currently assumes signal sampled at 1000 Hz !
%
% TODO: extend to signals not sampled at 1000 Hz
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, December 2017, updated November 2018
%-------------------------------------------------------------------------


%========== height threshold ==========
SigEpochs = (EEG(:,(PreStim+1):end)'>repmat(HeightThresh*std(EEG(:,1:PreStim-NotBaseline)'),PostStim+1,1))...
    +(EEG(:,(PreStim+1):end)'<repmat(-HeightThresh*std(EEG(:,1:PreStim-NotBaseline)'),PostStim+1,1));

%========== width threshold ==========
for jj = 1:size(SigEpochs,2)
    SigTimes = find(SigEpochs(:,jj));
    
    if ~isempty(SigTimes)
        x = diff(SigTimes')==1;
        f = find([false,x]~=[x,false]);
        d = f(2:2:end)-f(1:2:end-1);
        g = find(d>=WidthTresh);
        t = SigTimes(f(2*g-1));
        d2 = d(g);
        
        SigTimesBis = false(size(SigTimes));
        for jjj = 1:length(t)
            SigTimesBis(t(jjj):t(jjj)+d2(jjj))=true;
        end
        SigEpochs(:,jj) = dnif(SigTimesBis,PostStim+1);
    else
        SigEpochs(:,jj) = false(size(SigEpochs(:,jj)));
    end
end

Temp = EEG(:,(PreStim+1):end)-repmat(HeightThresh*std(EEG(:,1:PreStim-NotBaseline)'),PostStim+1,1)';
SigEpochsAmp = SigEpochs;
SigEpochsAmp(SigEpochsAmp==1)=Temp(SigEpochsAmp==1);

end