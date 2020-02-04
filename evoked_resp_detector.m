function [pklER,slER,rplER,msrlER,elER,dplER,msdlER,pkER,sER,rpER,msrER,eER,dpER,msdER,FinalMinPkPromB,FinalMinPkPromE] = evoked_resp_detector(EEG,SigMat,Forbidden)
% EVOKED_RESP_DETECTOR:  get peak and latency of average evoked response
%                 at peak, rising phase and onset. Responses are separated
%                 into "early" (0-50 ms) and "late" (50-500 ms) responses
%
% Usage:
%-------------------------------------------------------------------------
% [pklER,slER,rplER,msrlER,elER,dplER,msdlER,...
%   pkER,sER,rpER,msrER,eER,dpER,msdER,FinalMinPkPromB,...
%   FinalMinPkPromE] = evoked_resp_detector(EEG,SigSigmaMat,tScoreMat,Forbidden)
%
%
% Inputs:
%-------------------------------------------------------------------------
% EEG       : c x t EEG trace. Of importance: detecting response onset is
%             almost impossible at the trial level, but only very difficult
%             with a nice smoothed filtered average !
%             Therefore, variable EEG should an average or t-scores e.g.
%
% SigMat    : (t - b) x c significance binary matrix
%             (1: significant; 0: not significant)
%
% Forbidden : number of samples to ignore after stimulus onset before
%                searching for peak (even if these are significant in
%                SigMatStrict and SigMatLoose)
%
% c : number of channels
% t : number of time points
% b : number of pre-stimulus time points
%
%
% Outputs:
%-------------------------------------------------------------------------
% pkER            : ERP peak amplitude (c x n cell array)
%
% pklER           : ERP peak latency
%
% rpER            : amplitude of ERP at 50% rising phase
%
% rplER           : latency of ERP 50% rising phase
%
% msrER           : ERP maximal slope (between onset and peak)
%
% msrlER          : ERP maximal slope latency (between onset and peak)
%
% sER             : ERP amplitude at onset
%
% slER            : ERP onset latency
%
% eER             : ERP amplitude at offset
%
% elER            : ERP offset latency
%
% dpER            : amplitude of ERP at 50% descending phase
%
% dplER           : latency of ERP at 50% descending phase
%
% msdER           : ERP maximal slope (between peak and offset)
%
% msdlER          : ERP maximal slope (between peak and offset)
%
% FinalMinPkPromB : final minimal peak prominence for response onset
%
% FinalMinPkPromE : final minimal peak prominence for response offset
%
% n : number of peaks found
%
%
% Nota bene:
%-------------------------------------------------------------------------
% If you want to detect multiple peaks within an epoch,
% that is currently impossible with the current approach. That is
% the only drawback: this algorithm considers that each block is
% some sort of separate evoked response... Choose your significance
% level carefully!
%
% Currently assumes signal sampled at 1000 Hz !
%
% TODO: extend to signals not sampled at 1000 Hz (or just think in terms of
%       samples and beware with Initialization values below)
%
% RM@FBMlab: cumbersome function, should be refactored and more modular
% 
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, December 2017
%-------------------------------------------------------------------------

%% Init
Limit = 500; % zero-crossing should occur between stim onset and 500 ms post-stim
InitPkProm = 1; % if t-scores are provided and number of trials is around 20, 2 provides something close to significance threshold at p < 0.05
MPPdecreaseRate = 0.01; % rate at which PkProm should decrease if no peak is found

B = size(EEG,2)-size(SigMat,2);
if size(EEG,1)~=size(SigMat,1), error('Input size mismatch'), end;

% Exclude artifact period duration
SigMat(:,1:Forbidden)=0;

[pklER,slER,rplER,msrlER,elER,dplER,msdlER,pkER,sER,rpER,msrER,eER,dpER,msdER,FinalMinPkPromB,FinalMinPkPromE] = deal(cell(size(EEG,1),1)); % we do not know in advance how many peaks will be detected but we need to know if a channel elicited no response, even if it is in last channels !
for c = 1:size(EEG,1)
    if ~isempty(find(abs(diff(sign(EEG(c,B:B+Limit))))>1,1)) % try to detect peaks only if there is at least 1 zero-crossing between stimulation time and time Limit post-stim
        
        %% detect blocks of significance
        Blocks = diff(SigMat(c,:));
        Onsets = find(Blocks==1)+1;
        Offsets = find(Blocks==-1)+1;
        Nblocks1 = length(Onsets);
        Nblocks2 = length(Offsets);
        
        % sometimes Onsets or Offsets can be empty because only 1
        % event near the end or the beginning of the epoch...
        if isempty(Onsets)
            Onsets = 1;
        end
        if isempty(Offsets)
            Offsets = length(Blocks);
        end
        % sometimes there are multiple events and one of them is close to
        % the beginning or the end of the epoch...
        if Nblocks1 ~= Nblocks2
            if (Offsets(end)-Onsets(end))<1
                Offsets(end+1) = length(Blocks);
            elseif (Offsets(1)-Onsets(1))<1
                OnsetsBKP = Onsets;
                Onsets = nan(1,length(Onsets)+1);
                Onsets(1) = 1;
                Onsets(2:end) = OnsetsBKP;
            end
            Nblocks1 = length(Onsets);
            Nblocks2 = length(Offsets);
            if Nblocks1 ~= Nblocks2
                error('???')
            else
                Nblocks = Nblocks1;
            end
        else
            Nblocks = Nblocks1;
        end
        
        %% detect maximum / minimum within blocks
        % Nota bene: if you want to detect multiple peaks within an epoch,
        % that is currently impossible with the current approach. That is
        % the only drawback: this algorithm considers that each block is
        % some sort of separate evoked response... Choose your significance
        % level carefully!
        for b = 1:Nblocks
            
            if median(EEG(c,B+Onsets(b):B+Offsets(b)))<1
                [pkER{c,b},pklER{c,b}] = min(EEG(c,B+Onsets(b):B+Offsets(b)));
            elseif median(EEG(c,B+Onsets(b):B+Offsets(b)))>=1
                [pkER{c,b},pklER{c,b}] = max(EEG(c,B+Onsets(b):B+Offsets(b)));
            end
            
            pklER{c,b} = B+Onsets(b)+pklER{c,b}-1;
            pkER{c,b} = EEG(c,pklER{c,b});
            
        end
        
        IBBb = B;
        for b = 1:Nblocks
            IBBb(b+1)=Offsets(b)+B;
        end
        
        IBBe = Onsets+B;
        IBBe(end+1) = size(EEG,2);
        
        %========= BACKWARD DETECTION OF RESPONSE ONSET, RISING PHASE AND MAX SLOPE =========
        for b = 1:Nblocks
            
            PkProm = InitPkProm;
            FinalMinPkPromB{c,b} = PkProm;
            
            %% find first deflection before peak or zero-crossing in between
            % detection of response onset is crucial as it will also define
            % the rising phase and the maximal slope => this algorithm
            % considers that response onset should start before
            % significance, i.e. that evoked response's onset is
            % subthreshold.
            if length(IBBb(b):IBBe(b))<3
                slER{c,b} = IBBb(b);
                sER{c,b} = EEG(c,slER{c,b});
            else
                
                % before trying to find zero-crossing, get first deflection
                % of "reasonable" prominence, i.e. = 1
                if pkER{c,b}<1
                    [sER{c,b},slER{c,b}] = findpeaks(fliplr(EEG(c,IBBb(b):IBBe(b))),'npeaks',1,'minpeakprominence',PkProm);
                elseif pkER{c,b}>=1
                    [sER{c,b},slER{c,b}] = findpeaks(fliplr(-EEG(c,IBBb(b):IBBe(b))),'npeaks',1,'minpeakprominence',PkProm);
                end
                
                if ~isempty(slER{c,b})
                    %   ZC = find(abs(diff(sign(EEG(c,IBBb(b):IBBe(b)))))>1);
                    ZC = find(abs(diff(sign(EEG(c,IBBb(b):(IBBe(b)-slER{c,b}+1)))))>1);
                    
                    if ~isempty(ZC)
                        slER{c,b} = IBBb(b)+ZC(end);
                        sER{c,b} = EEG(c,slER{c,b});
                    else
                        slER{c,b} = IBBe(b)-slER{c,b}+1;
                        sER{c,b} = EEG(c,slER{c,b});
                    end
                    
                else
                    % in this case we lower minpeakprominence until we find
                    % even the smallest deflection
                    if length(IBBb(b):IBBe(b))<3
                        slER{c,b} = 1;
                        sER{c,b} = EEG(c,IBBb(b));
                    else
                        slER{c,b} = [];
                        while PkProm>=(MPPdecreaseRate) && isempty(slER{c,b})
                            
                            PkProm = PkProm-MPPdecreaseRate;
                            if pkER{c,b}<1
                                [sER{c,b},slER{c,b}] = findpeaks(fliplr(EEG(c,IBBb(b):IBBe(b))),'npeaks',1,'minpeakprominence',PkProm);
                            elseif pkER{c,b}>=1
                                [sER{c,b},slER{c,b}] = findpeaks(fliplr(-EEG(c,IBBb(b):IBBe(b))),'npeaks',1,'minpeakprominence',PkProm);
                            end
                            
                        end
                    end
                    
                    FinalMinPkProm{c,b} = PkProm;
                    
                    if isempty(slER{c,b}) % still empty but reached minpeakprominence of 0
                        
                        % AFAIK, this situation occurs only if the
                        % beginning of the epoch (i.e. stimulation onset)
                        % but the response started slightly earlier...
                        
                        slER{c,b} = B;
                        sER{c,b} = EEG(c,slER{c,b});
                        
                    else
                        
                        slER{c,b} = IBBe(b)-slER{c,b}+1;
                        sER{c,b} = EEG(c,slER{c,b});
                        
                    end
                end
                %                 if pkER{c,b}<1
                %                     [sER{c,b},slER{c,b}] = findpeaks(fliplr(EEG(c,IBBb(b):IBBe(b))),'npeaks',1,'minpeakprominence',PkProm);
                %                 elseif pkER{c,b}>=1
                %                     [sER{c,b},slER{c,b}] = findpeaks(fliplr(-EEG(c,IBBb(b):IBBe(b))),'npeaks',1,'minpeakprominence',PkProm);
                %                 end
            end
            
            %             if ~isempty(slER{c,b})
            %                 slER{c,b} = IBBe(b)-slER{c,b}+1;
            %                 sER{c,b} = EEG(c,slER{c,b});
            %
            %                 % if zero-crossing in between, use it as onset
            %                 ZC = find(abs(diff(sign(EEG(c,slER{c,b}:pklER{c,b}))))>1);
            %                 if ~isempty(ZC)
            %                     slER{c,b} = slER{c,b}+ZC(1);
            %                     sER{c,b} = EEG(c,slER{c,b});
            %                 end
            %
            %             else
            % in some cases the signal does not stabilize and just
            % decreases more or less continuously => this happens
            % especially between two big deflections of oppostite
            % polarity => in this case the zero-crossing points will
            % provide an accurate onset for the next response
            
            %             end
            
            %% now get 50% rising phase
            rplER{c,b} = round(mean([slER{c,b},pklER{c,b}]));
            rpER{c,b} = EEG(c,rplER{c,b});
            
            %% and max slope
            if median(EEG(c,slER{c,b}:pklER{c,b}))<1
                [msrER{c,b}, msrlER{c,b}] = min(diff(EEG(c,slER{c,b}:pklER{c,b})));
            elseif median(EEG(c,slER{c,b}:pklER{c,b}))>=1
                [msrER{c,b}, msrlER{c,b}] = max(diff(EEG(c,slER{c,b}:pklER{c,b})));
            end
            msrlER{c,b} = msrlER{c,b} + slER{c,b};
            msrER{c,b} = EEG(c,msrlER{c,b});
            
        end
        
        %========= FORWARD DETECTION OF RESPONSE OFFSET, DESCENDING PHASE AND MAX SLOPE =========
        for b = 1:Nblocks
            
            PkProm = InitPkProm;
            FinalMinPkPromE{c,b} = PkProm;
            
            %% find first deflection after peak (same rationale as for response onset)
            if length(IBBb(b+1):IBBe(b+1))<3
                elER{c,b} = IBBe(b+1);
                eER{c,b} = EEG(c,elER{c,b});
            else
                
                % before trying to find zero-crossing, get first deflection
                % of "reasonable" prominence, i.e. = 1
                if pkER{c,b}<1
                    [eER{c,b},elER{c,b}] = findpeaks(EEG(c,IBBb(b+1):IBBe(b+1)),'npeaks',1,'minpeakprominence',PkProm);
                elseif pkER{c,b}>=1
                    [eER{c,b},elER{c,b}] = findpeaks(-EEG(c,IBBb(b+1):IBBe(b+1)),'npeaks',1,'minpeakprominence',PkProm);
                end
                
                if ~isempty(elER{c,b})
                    %   ZC = find(abs(diff(sign(EEG(c,IBBb(b+1):IBBe(b+1)))))>1);
                    ZC = find(abs(diff(sign(EEG(c,IBBb(b+1):(IBBb(b+1)+elER{c,b}+1)))))>1);
                    if ~isempty(ZC)
                        elER{c,b} = IBBb(b+1)+ZC(1);
                        eER{c,b} = EEG(c,elER{c,b});
                    else
                        elER{c,b} = IBBb(b+1)+elER{c,b};
                        eER{c,b} = EEG(c,elER{c,b});
                    end
                    
                else
                    % in this case we lower minpeakprominence until we find
                    % even the smallest deflection
                    if length(IBBb(b+1):IBBe(b+1))<3
                        elER{c,b} = IBBe(b+1);
                        eER{c,b} = EEG(c,elER{c,b});
                    else
                        elER{c,b} = [];
                        while PkProm>=(MPPdecreaseRate) && isempty(elER{c,b})
                            
                            PkProm = PkProm-MPPdecreaseRate;
                            if pkER{c,b}<1
                                [eER{c,b},elER{c,b}] = findpeaks(EEG(c,IBBb(b+1):IBBe(b+1)),'npeaks',1,'minpeakprominence',PkProm);
                            elseif pkER{c,b}>=1
                                [eER{c,b},elER{c,b}] = findpeaks(-EEG(c,IBBb(b+1):IBBe(b+1)),'npeaks',1,'minpeakprominence',PkProm);
                            end
                            
                        end
                    end
                    
                    FinalMinPkProm{c,b} = PkProm;
                    
                    if isempty(elER{c,b}) % still empty but reached minpeakprominence of 0
                        
                        % AFAIK, it means that we reached end of epoch
                        
                        elER{c,b} = IBBe(end);
                        eER{c,b} = EEG(c,elER{c,b});
                        
                    else
                        
                        elER{c,b} = IBBb(b+1)+elER{c,b};
                        eER{c,b} = EEG(c,elER{c,b});
                        
                    end
                end
                %                 if pkER{c,b}<1
                %                     [eER{c,b},elER{c,b}] = findpeaks(EEG(c,IBBb(b+1):IBBe(b+1)),'npeaks',1,'minpeakprominence',PkProm);
                %                 elseif pkER{c,b}>=1
                %                     [eER{c,b},elER{c,b}] = findpeaks(-EEG(c,IBBb(b+1):IBBe(b+1)),'npeaks',1,'minpeakprominence',PkProm);
                %                 end
            end
            
            %             if ~isempty(elER{c,b})
            %                 elER{c,b} = IBBb(b+1)+elER{c,b}+1;
            %                 eER{c,b} = EEG(c,elER{c,b});
            %
            %                 % if zero-crossing in between, use it as offset
            %                 ZC = find(abs(diff(sign(EEG(c,pklER{c,b}:elER{c,b}))))>1);
            %                 if ~isempty(ZC)
            %                     elER{c,b} = pklER{c,b}+ZC(1);
            %                     eER{c,b} = EEG(c,elER{c,b});
            %                 end
            %
            %             else
            
            %             end
            
            %% now get 50% descending phase
            dplER{c,b} = round(mean([elER{c,b},pklER{c,b}]));
            dpER{c,b} = EEG(c,dplER{c,b});
            
            %% and max slope
            if median(EEG(c,pklER{c,b}:elER{c,b}))<1
                [msdER{c,b}, msdlER{c,b}] = max(diff(EEG(c,pklER{c,b}:elER{c,b})));
            elseif median(EEG(c,pklER{c,b}:elER{c,b}))>=1
                [msdER{c,b}, msdlER{c,b}] = min(diff(EEG(c,pklER{c,b}:elER{c,b})));
            end
            msdlER{c,b} = msdlER{c,b} + pklER{c,b};
            msdER{c,b} = EEG(c,msdlER{c,b});
            
        end
    end
end

