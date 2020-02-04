function [pkER,pklER,pkwER,pkpER,msER,mslER,sER,slER,pkLR,pklLR,pkwLR,pkpLR,msLR,mslLR,sLR,slLR] = get_CCEP_resp(avgEEG,SigMatStrict,SigMatLoose,Forbidden)
% GET_CCEP_RESP:  get peak and latency of average evoked response
%                 at peak, rising phase and onset. Responses are separated
%                 into "early" (0-50 ms) and "late" (50-500 ms) responses
%
% Usage:
%-------------------------------------------------------------------------
% [pkER,pklER,pkwER,pkpER,msER,mslER,sER,slER,pkLR,pklLR,pkwLR,...
%   pkpLR,msLR,mslLR,sLR,...
%   slLR] = get_CCEP_resp(avgEEG,SigMatStrict,SigMatLoose,Forbidden)
%
%
% Inputs:
%-------------------------------------------------------------------------
% avgEEG       : c x t average evoked response
%
% SigMatStrict : (t - b) x c significance binary matrix
%                (1: significant; 0: not significant), with conservative
%                significance threshold
%
% SigMatLoose  : (t - b) x c significance binary matrix
%                (1: significant; 0: not significant), with liberal
%                significance threshold
%
% Forbidden    : number of samples to ignore after stimulus onset before
%                searching for peak (even if these are significant in
%                SigMatStrict and SigMatLoose)
%
% c = number of channels
% t = number of time points
% b = number of pre-stimulus time points
%
%
% Outputs:
%-------------------------------------------------------------------------
% pkER  : early response peak amplitude
%
% pklER : early response peak latency
%
% pkwER : early response peak width at half prominence
%
% pkpER : early response peak prominence
%
% msER  : early response maximal slope (rising phase)
%
% mslER : early response maximal slope latency
%
% sER   : early response amplitude at start
%
% slER  : early response start latency
%
% pkLR  : late response peak amplitude
%
% pklLR : late response peak latency
%
% pkwLR : late response peak width at half prominence
%
% pkpLR : late response peak prominence
%
% msLR  : late response maximal slope (rising phase)
%
% mslLR : late response maximal slope latency
%
% sLR   : late response amplitude at start
%
% slLR  : late response start latency
%
% Nota bene:
%-------------------------------------------------------------------------
% Currently assumes signal sampled at 1000 Hz !
%
% TODO: extend to signals not sampled at 1000 Hz
%
% RM@FBMlab: cumbersome function, should be refactored and more modular
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, November 2017
%-------------------------------------------------------------------------

% Init
Early = 51;
Late = 501;
Limit = 250; % zero-crossing should occur between stim onset and 250 ms post-stim

B = size(avgEEG,2)-size(SigMatStrict,1);
check_input_args(avgEEG,SigMatStrict,SigMatLoose,Forbidden,Early);

pkER = nan(size(avgEEG,1),1);
pklER = nan(size(avgEEG,1),1);
pkwER = nan(size(avgEEG,1),1);
pkpER = nan(size(avgEEG,1),1);
msER = nan(size(avgEEG,1),1);
mslER = nan(size(avgEEG,1),1);
sER = nan(size(avgEEG,1),1);
slER = nan(size(avgEEG,1),1);

pkLR = nan(size(avgEEG,1),1);
pklLR = nan(size(avgEEG,1),1);
pkwLR = nan(size(avgEEG,1),1);
pkpLR = nan(size(avgEEG,1),1);
msLR = nan(size(avgEEG,1),1);
mslLR = nan(size(avgEEG,1),1);
sLR = nan(size(avgEEG,1),1);
slLR = nan(size(avgEEG,1),1);

% find zero-crossings function:
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);

for c = 1:size(avgEEG,1)
    if ~isempty(zci(avgEEG(c,B:B+Limit))) % try to detect peaks only if there is at least 1 zero-crossing between stimulation time and 500 post-stim
        if any(SigMatStrict((Forbidden+1):Early,c)) % early response?
            % detect response peak
            if median(avgEEG(c,B+find(SigMatStrict((Forbidden+1):Early,c))))>=0
                temp = findpeaks(avgEEG(c,B+((Forbidden+1):Early))','npeaks',1,'minpeakprominence',3);
                if ~isempty(temp)
                    [pkER(c),pklER(c),pkwER(c),pkpER(c)] = findpeaks(avgEEG(c,B+((Forbidden+1):Early))','npeaks',1,'minpeakprominence',3);
                else
                    pkER(c) = nan;
                    pklER(c) = nan;
                    pkwER(c) = nan;
                    pkpER(c) = nan;
                end
            else
                temp = findpeaks(-avgEEG(c,B+((Forbidden+1):Early))','npeaks',1,'minpeakprominence',3);
                if ~isempty(temp)
                    [pkER(c),pklER(c),pkwER(c),pkpER(c)] = findpeaks(-avgEEG(c,B+((Forbidden+1):Early))','npeaks',1,'minpeakprominence',3);
                else
                    pkER(c) = nan;
                    pklER(c) = nan;
                    pkwER(c) = nan;
                    pkpER(c) = nan;
                end
            end
            if ~isnan(pkER(c))
                
                % detect response max slope:
                FirstSigBefore = find(abs(diff(SigMatLoose((pklER(c)):-1:1,c))),1);
                if isempty(FirstSigBefore) || FirstSigBefore < (B+1)
                    FirstSigBefore = (B+1);
                end
                % could use also half prominence but
                % probably better to look at max slope...
                pklER(c) = (B+Forbidden)+pklER(c);
                pkER(c) = avgEEG(c,pklER(c));
                if median(avgEEG(c,B+find(SigMatStrict((Forbidden+1):Early,c))))>=0
                    [msER(c),mslER(c)] = max(diff(avgEEG(c,FirstSigBefore:pklER(c))));
                else
                    [msER(c),mslER(c)] = min(diff(avgEEG(c,FirstSigBefore:pklER(c))));
                end
                
                mslER(c) = FirstSigBefore+mslER(c);
                msER(c) = avgEEG(c,mslER(c));
                
                % detect valley (response start) after first significant amplitude
                % and before max slope:
                if median(avgEEG(c,B+find(SigMatStrict((Forbidden+1):Early,c))))>=0
                    if abs(mslER(c) - FirstSigBefore) < 2
                        sER(c) = nan;
                        slER(c) = nan;
                    else
                        if ~isempty(findpeaks(-(avgEEG(c,FirstSigBefore:mslER(c))),'npeaks',1,'minpeakprominence',3))
                            [sER(c),slER(c)] = findpeaks(-(avgEEG(c,FirstSigBefore:mslER(c))),'npeaks',1,'minpeakprominence',3);
                            slER(c) = slER(c)+FirstSigBefore;
                            sER(c) = avgEEG(c,slER(c));
                        else
                            slER(c) = FirstSigBefore;
                            sER(c) = avgEEG(c,FirstSigBefore);
                        end
                    end
                else
                    if abs(mslER(c) - FirstSigBefore) < 2
                        sER(c) = nan;
                        slER(c) = nan;
                    else
                        if ~isempty(findpeaks((avgEEG(c,FirstSigBefore:mslER(c))),'npeaks',1,'minpeakprominence',3))
                            [sER(c),slER(c)] = findpeaks((avgEEG(c,FirstSigBefore:mslER(c))),'npeaks',1,'minpeakprominence',3);
                            slER(c) = slER(c)+FirstSigBefore;
                            sER(c) = avgEEG(c,slER(c));
                        else
                            slER(c) = FirstSigBefore;
                            sER(c) = avgEEG(c,FirstSigBefore);
                        end
                    end
                end
            else
                msER(c) = nan;
                mslER(c) = nan;
                sER(c) = nan;
                slER(c) = nan;
            end
        end
        if any(SigMatStrict(Early:Late,c)) % late response?
            % detect response peak
            if median(avgEEG(c,B+find(SigMatStrict(Early:Late,c))))>=0
                temp = findpeaks(avgEEG(c,B+(Early:Late))','npeaks',1,'minpeakprominence',3);
                if ~isempty(temp)
                    [pkLR(c),pklLR(c),pkwLR(c),pkpLR(c)] = findpeaks(avgEEG(c,B+(Early:Late))','npeaks',1,'minpeakprominence',3);
                else
                    pkLR(c) = nan;
                    pklLR(c) = nan;
                    pkwLR(c) = nan;
                    pkpLR(c) = nan;
                end
            else
                temp = findpeaks(-avgEEG(c,B+(Early:Late))','npeaks',1,'minpeakprominence',3);
                if ~isempty(temp)
                    [pkLR(c),pklLR(c),pkwLR(c),pkpLR(c)] = findpeaks(-avgEEG(c,B+(Early:Late))','npeaks',1,'minpeakprominence',3);
                else
                    pkLR(c) = nan;
                    pklLR(c) = nan;
                    pkwLR(c) = nan;
                    pkpLR(c) = nan;
                end
            end
            if ~isnan(pkLR(c))
                
                % detect response max slope:
                FirstSigBefore = find(abs(diff(SigMatLoose((pklLR(c)):-1:1,c))),1);
                if isempty(FirstSigBefore) || FirstSigBefore < (B+Early)
                    FirstSigBefore = (B+Early);
                end
                % could use also half prominence but
                % probably better to look at max slope...
                pklLR(c) = (B+Early-1)+pklLR(c);
                pkLR(c) = avgEEG(c,pklLR(c));
                if median(avgEEG(c,B+find(SigMatStrict(Early:Late,c))))>=0
                    [msLR(c),mslLR(c)] = max(abs(diff(avgEEG(c,FirstSigBefore:pklLR(c)))));
                else
                    [msLR(c),mslLR(c)] = min(abs(diff(avgEEG(c,FirstSigBefore:pklLR(c)))));
                end
                
                mslLR(c) = FirstSigBefore+mslLR(c);
                msLR(c) = avgEEG(c,mslLR(c));
                
                % detect valley (response start) after first significant amplitude
                % and before max slope:
                if median(avgEEG(c,B+find(SigMatStrict(Early:Late,c))))>=0
                    if abs(mslLR(c) - FirstSigBefore) < 2
                        sLR(c) = nan;
                        slLR(c) = nan;
                    else
                        if ~isempty(findpeaks(-(avgEEG(c,FirstSigBefore:mslLR(c))),'npeaks',1,'minpeakprominence',3))
                            [sLR(c),slLR(c)] = findpeaks(-(avgEEG(c,FirstSigBefore:mslLR(c))),'npeaks',1,'minpeakprominence',3);
                            slLR(c) = slLR(c)+FirstSigBefore;
                            sLR(c) = avgEEG(c,slLR(c));
                        else
                            slLR(c) = FirstSigBefore;
                            sLR(c) = avgEEG(c,FirstSigBefore);
                        end
                    end
                else
                    if abs(mslLR(c) - FirstSigBefore) < 2
                        sLR(c) = nan;
                        slLR(c) = nan;
                    else
                        if ~isempty(findpeaks((avgEEG(c,FirstSigBefore:mslLR(c))),'npeaks',1,'minpeakprominence',3))
                            [sLR(c),slLR(c)] = findpeaks((avgEEG(c,FirstSigBefore:mslLR(c))),'npeaks',1,'minpeakprominence',3);
                            slLR(c) = slLR(c)+FirstSigBefore;
                            sLR(c) = avgEEG(c,slLR(c));
                        else
                            slLR(c) = FirstSigBefore;
                            sLR(c) = avgEEG(c,FirstSigBefore);
                        end
                    end
                end
            else
                msLR(c) = nan;
                mslLR(c) = nan;
                sLR(c) = nan;
                slLR(c) = nan;
            end
        end
    end
end

end

%########################## INTERNAL FUNCTIONS ###########################
function check_input_args(avgEEG,SigMatStrict,SigMatLoose,Forbidden,Early)

test(1) = size(SigMatStrict,1)~=size(SigMatLoose,1);
test(2) = size(SigMatStrict,2)~=size(SigMatLoose,2);
test(3) = size(avgEEG,1)~=size(SigMatStrict,2);
test(4) = size(avgEEG,1)~=size(SigMatLoose,2);

if sum(test)>0
    error('Input dimensions mismatch');
end

if Forbidden>=(Early-1)
    error('Forbidden period is too long to analyze early responses')
end
end
