function [gfpHRP,gfpHRPval,GFP,gfpPeak,gfpPeakVal,WinSize] = find_GFP_HRP(EEG,Origin,WinSize)
% find_GFP_HRP: find time point of 50% rising phase
%
% [gfpHRP,gfpHRPval,GFP,gfpPeak,gfpPeakVal] = find_GFP_HRP(EEG,Origin,WinSize)
%
%  Inputs
% --------
% EEG: channels x time EEG traces
% Origin: time point of stimulus
%         (cannot be < 1 but is usually bigger than 1 because of pre-stimulus period)
% WinSize: window size around Origin for peak detection
%
%  Outputs
% ---------
% gfpHRP: time point of 50% rising phase
% gfpHRPval: GFP value at 50% rising phase
% GFP: global field power
% gfpPeak: time point of peak
% gfpPeakVal: GFP value at peak
% WinSize: (might be updated if algorithm struggles to find the peak...)
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, April 2018, last updated August 2018
%-------------------------------------------------------------------------

if Origin<1
    error('Origin cannot be less than 1 time frame (do not forget to consider also pre-stimulus period)')
end
if nargin<3
    WinSize = 200;
end
% OrigWinSize = WinSize;

GFP = computegfp(EEG');

% [MaxGFP,gfpPeak] = max(GFP(Origin-round(WinSize/2):Origin+round(WinSize/2)));
[MaxGFP,gfpPeak] = max(GFP(Origin:Origin+round(WinSize)));
% [MaxGFP,gfpPeak] = max(GFP(Origin:Origin+round(WinSize)));
% gfpPeak = Origin-round(WinSize/2)+gfpPeak-1;
gfpPeak = Origin+gfpPeak-1;
% gfpPeak = Origin+gfpPeak-1;
BKPpeak = gfpPeak;
FlagDoNotShiftAgain = false;
if ~(GFP(gfpPeak)>=GFP(gfpPeak+1) && GFP(gfpPeak)>=GFP(gfpPeak-1))
    warning('Peak might be on the edge, trying findpeaks instead...')
%     [MaxGFP,gfpPeak] = findpeaks(GFP(Origin-round(WinSize/2):Origin+round(WinSize/2)));
    [MaxGFP,gfpPeak] = findpeaks(GFP(Origin:Origin+round(WinSize)));
    
    if isempty(gfpPeak)
        FlagDoNotShiftAgain = true;
        gfpPeak = BKPpeak;
        while ~(GFP(gfpPeak)>=GFP(gfpPeak+1) && GFP(gfpPeak)>=GFP(gfpPeak-1))
            WinSize = WinSize+1;
%             [MaxGFP,gfpPeak] = max(GFP(Origin-round(WinSize/2):Origin+round(WinSize/2)));
            [MaxGFP,gfpPeak] = max(GFP(Origin:Origin+round(WinSize)));
%             gfpPeak = Origin-round(WinSize/2)+gfpPeak-1;
            gfpPeak = Origin+gfpPeak-1;
        end
    end
    
    % there might be multiple peaks...
    if length(gfpPeak)>1
        clear DistO
        for i = 1:length(gfpPeak)
            DistO(i) = abs(gfpPeak(i)-Origin);
        end
        [~,idx] = min(DistO);
        gfpPeak = gfpPeak(idx);
        MaxGFP = MaxGFP(idx);
    end
    if ~FlagDoNotShiftAgain
%         gfpPeak = Origin-round(WinSize/2)+gfpPeak-1;
        gfpPeak = Origin+gfpPeak-1;
    end
end

% estimating Half-Rising Phase (HRP)
BelowHRP = gfpPeak-find(GFP(gfpPeak:-1:1)<(MaxGFP/2),1)+1;
if isempty(BelowHRP)
    BelowHRP = 1;
end

% if BelowHRP<(Origin-round(OrigWinSize/2))
if BelowHRP<=(Origin-round(WinSize)) % should not happen, except if response amplitude is not so high...
%     % this can happen when the amplitude of the response is not so high...
%     warning('GFP amplitude at peak looks not high enough to decrease enough soon enough... Accepting 45% rising phase...')
%     BelowHRP = gfpPeak-find(GFP(gfpPeak:-1:1)<(MaxGFP/100*45),1)+1;
    MedianGFPvalBefore = median(GFP(gfpPeak:-1:1));
    MidGFP = MedianGFPvalBefore + ((MaxGFP - MedianGFPvalBefore)/2);
    BelowHRP = gfpPeak-find(GFP(gfpPeak:-1:1)<MidGFP,1)+1;
end
% if BelowHRP<(Origin-round(OrigWinSize/2))
%     warning('GFP amplitude looks really low... Using findpeaks... Check GFP % at HRP!')
%     [~,BelowHRPtemp] = findpeaks(-GFP(gfpPeak:-1:1));
%     if ~isempty(BelowHRPtemp)
%     BelowHRP = gfpPeak-BelowHRPtemp(1);
%     else
%         error('Nor GFP 50%, neither GFP valley before peak were found! Check your data!')
%     end
% end

AboveHRP = BelowHRP+1;
[~,Closest2HRP] = min([abs((MaxGFP/2)-GFP(BelowHRP)),abs((MaxGFP/2)-GFP(AboveHRP))]);
asdf = [BelowHRP,AboveHRP];
gfpHRP = asdf(Closest2HRP);

% Time = find(diff(find(sign(single((GFP(Origin:-1:1)-GFP(Origin)/2)>0))))>1);
% Time = Origin-Time(1)+1;
gfpHRPval = GFP(gfpHRP);
gfpPeakVal = GFP(gfpPeak);

end


