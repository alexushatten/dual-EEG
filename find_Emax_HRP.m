function [HRP,HRPval,Emax,PeakEmax,PeakVal,WinSize] = find_Emax_HRP(EEG,Origin,WinSize)
% find_Emax_HRP: find time point of 50% rising phase
%
% [HRP,HRPval,Emax,Peak,PeakVal] = find_Emax_HRP(EEG,Origin,WinSize)
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
% HRP: time point of 50% rising phase
% HRPval: Emax value at 50% rising phase
% Emax: electrode with highest (absolute) amplitude at peak
% Peak: time point of peak
% PeakVal: Emax value at peak
% WinSize: (might be updated if algorithm struggles to find the peak...)
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, May 2018, last updated August 2018
%-------------------------------------------------------------------------

if Origin<1
    error('Origin cannot be less than 1 time frame (do not forget to consider also pre-stimulus period)')
end
if nargin<3
    WinSize = 200;
end
OrigWinSize = WinSize;

GFP = computegfp(EEG');

% [~,Peak] = max(GFP(Origin-round(WinSize/2):Origin+round(WinSize/2)));
[~,Peak] = max(GFP(Origin:Origin+round(WinSize)));
% [~,Peak] = max(GFP(Origin:Origin+round(WinSize)));
% Peak = Origin-round(WinSize/2)+Peak-1;
Peak = Origin+Peak-1;
% Peak = Origin+Peak-1;
BKPpeak = Peak;
FlagDoNotShiftAgain = false;
if ~(GFP(Peak)>=GFP(Peak+1) && GFP(Peak)>=GFP(Peak-1))
    warning('Preliminary GFP peak might be on the edge, trying findpeaks instead...')
%     [~,Peak] = findpeaks(GFP(Origin-round(WinSize/2):Origin+round(WinSize/2)));
    [~,Peak] = findpeaks(GFP(Origin:Origin+round(WinSize)));
    
    if isempty(Peak)
        FlagDoNotShiftAgain = true;
        Peak = BKPpeak;
        while ~(GFP(Peak)>=GFP(Peak+1) && GFP(Peak)>=GFP(Peak-1))
            WinSize = WinSize+1;
%             [~,Peak] = max(GFP(Origin-round(WinSize/2):Origin+round(WinSize/2)));
            [~,Peak] = max(GFP(Origin:Origin+round(WinSize)));
%             Peak = Origin-round(WinSize/2)+Peak-1;
            Peak = Origin+Peak-1;
        end
    end
    
    % there might be multiple peaks...
    if length(Peak)>1
        clear DistO
        for i = 1:length(Peak)
            DistO(i) = abs(Peak(i)-Origin);
        end
        [~,idx] = min(DistO);
        Peak = Peak(idx);
    end
    if ~FlagDoNotShiftAgain
%         Peak = Origin-round(WinSize/2)+Peak-1;
        Peak = Origin+Peak-1;
    end
end

[~,Emax] = max(abs(EEG(:,Peak)));
% in next line we further restrain the possible options because sometimes
% it could lead to a jump that is too big... in any case peak of Emax
% should normally be close to GFP max...
if EEG(Emax,Peak)>0
    [MaxEmax,PeakEmax] = max(EEG(Emax,Peak-round(OrigWinSize/2):Peak+round(OrigWinSize/2)));
elseif EEG(Emax,Peak)<0
    [MaxEmax,PeakEmax] = min(EEG(Emax,Peak-round(OrigWinSize/2):Peak+round(OrigWinSize/2)));
end
% [MaxEmax,PeakEmax] = max(abs(EEG(Emax,Peak:Peak+round(OrigWinSize))));
PeakEmax = Peak-round(OrigWinSize/2)+PeakEmax-1;
MaxEmax = EEG(Emax,PeakEmax);
% PeakEmax = Peak+PeakEmax-1;
BKPpeak = PeakEmax;
FlagDoNotShiftAgain = false;
if ~(abs(EEG(Emax,PeakEmax))>=abs(EEG(Emax,PeakEmax+1)) && abs(EEG(Emax,PeakEmax))>=abs(EEG(Emax,PeakEmax-1)))
    warning('Peak might be on the edge, check results!')
    [MaxEmax,PeakEmax] = findpeaks(abs(EEG(Emax,Peak-round(OrigWinSize/2):Peak+round(OrigWinSize/2))));
%     [MaxEmax,PeakEmax] = findpeaks(abs(EEG(Emax,Peak:Peak+round(OrigWinSize))));
    MaxEmax = EEG(Emax,PeakEmax);

    if isempty(PeakEmax)
        FlagDoNotShiftAgain = true;
        PeakEmax = BKPpeak;
        WinSize = OrigWinSize;
        while ~(abs(EEG(PeakEmax))>=abs(EEG(PeakEmax+1)) && abs(EEG(PeakEmax))>=abs(EEG(PeakEmax-1)))
            WinSize = WinSize+1;
%             [~,PeakEmax] =
%             max(abs(EEG(Origin-round(WinSize/2):Origin+round(WinSize/2))));       % => THIS WAS WRONG!!!
            [MaxEmax,PeakEmax] = max(abs(EEG(Peak-round(OrigWinSize/2):Peak+round(WinSize/2))));
%             [~,PeakEmax] = max(abs(EEG(Peak:Peak+round(WinSize))));
%             PeakEmax = Origin-round(WinSize/2)+PeakEmax-1;       % => THIS WAS WRONG!!!
            PeakEmax = Peak-round(OrigWinSize/2)+PeakEmax-1;
%             PeakEmax = Peak+PeakEmax-1;
        end
    end
    
    % there might be multiple peaks...
    if length(PeakEmax)>1
        clear DistO
        for i = 1:length(PeakEmax)
            DistO(i) = abs(PeakEmax(i)-Origin);
        end
        [~,idx] = min(DistO); % we keep the one closest to the origin! (should not be negative anyhow)
        if min(DistO)<0 % just to be sure
            error('Peak is before Origin!')
        end
        PeakEmax = PeakEmax(idx);
    end
    if ~FlagDoNotShiftAgain
        PeakEmax = Peak-round(OrigWinSize/2)+PeakEmax-1;
%         PeakEmax = Peak+PeakEmax-1;
    end
end

if numel(PeakEmax)>1 % this should not happen anyhow, just chckin'
%     warning('Multiple peaks found!... Check your output!')
    error('Multiple peaks found!... Check your code!!')
end

if MaxEmax>0
    BelowHRP = PeakEmax-find(EEG(Emax,PeakEmax:-1:1)<(MaxEmax/2),1)+1;
elseif MaxEmax<0
    BelowHRP = PeakEmax-find(EEG(Emax,PeakEmax:-1:1)>(MaxEmax/2),1)+1;
else
    error('Inconsistent value of maximum for electrode with maximal response')
end
if isempty(BelowHRP)
    BelowHRP = 1;
end
% if BelowHRP<(Origin-round(OrigWinSize/2))
if BelowHRP<=(Origin-round(WinSize))
%     % this can happen when the amplitude of the response is not so high...
%     warning('Emax amplitude at peak looks not high enough to decrease enough soon enough... Accepting 45% rising phase...')
%     BelowHRP = PeakEmax-find(abs(EEG(Emax,PeakEmax:-1:1))<(MaxEmax/100*45),1)+1;
    MedianEEGvalBefore = median(EEG(Emax,PeakEmax:-1:1));
    MidEEG = MedianEEGvalBefore + ((MaxEmax - MedianEEGvalBefore)/2);
    if MidEEG>0
        BelowHRP = PeakEmax-find(EEG(Emax,PeakEmax:-1:1)<MidEEG,1)+1;
    elseif MidEEG<0
        BelowHRP = PeakEmax-find(EEG(Emax,PeakEmax:-1:1)>MidEEG,1)+1;
    end
end
% if BelowHRP<(Origin-round(OrigWinSize/2))
%     warning('Emax amplitude looks really low... Using findpeaks... Check Emax % at HRP!')
%     [~,BelowHRPtemp] = findpeaks(-abs(EEG(Emax,PeakEmax:-1:1)));
%     if ~isempty(BelowHRPtemp)
%     BelowHRP = PeakEmax-BelowHRPtemp(1);
%     else
%         error('Nor Emax half-rising phase, neither Emax valley before peak were found! Check your data!')
%     end
% end

AboveHRP = BelowHRP+1;
[~,Closest2HRP] = min([abs((MaxEmax/2)-abs(EEG(Emax,BelowHRP))),abs((MaxEmax/2)-abs(EEG(Emax,AboveHRP)))]);
asdf = [BelowHRP,AboveHRP];
HRP = asdf(Closest2HRP);

% Time = find(diff(find(sign(single((GFP(Origin:-1:1)-GFP(Origin)/2)>0))))>1);
% Time = Origin-Time(1)+1;
HRPval = EEG(Emax,HRP);
PeakVal = EEG(Emax,PeakEmax);

end


