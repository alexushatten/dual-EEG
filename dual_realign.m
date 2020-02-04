function dual_realign(filename1,filename2, mkrfile, fs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GOAL: Realign the extracranial and intracranial recordings after resampling
% and interpolating data at 1000 Hz
%
% NOTE:
% - Assumes resampling frequency matches hdEEG sampling frequency, i.e. 1000Hz
% - Prerequisite: need exact same number and matching marker labels (tested
% first), use mkr = dual_clean_mkr
% - The sampling frequency is empirically based on the median period for
% each recording, e.g. if set fs is 2048 Hz, it may end up being 2049 Hz if
% most frames encompass 2049 datapoints.
%
% CALLING FUNCTIONS
% - [data des] = dual_load_sef(filname)
% - mkr = dual_clean_mkr(textfilenames)
%
% INPUT in prompt:
% - INPUT DIRECTORY contains file 1 (icEEG), file 2 (hdEEG) and merged marker file
% - filename1, data 1:     intracranial data (micromed 2048Hz)
% - filename2, data 2:     scalp data (EGI 1000 Hz)
% - marker file:    markers with labels (col 1), and datapoints for icEEG
%                   (col 2) and hdEEG (col 3), obtained from dual_clean_mrk
% - start/stop: optional, default 1:end, markers label at which to start
%               and stop
% - fs,resampf:    optional, resampling frequency, 1000 Hz by default
%
% OUTPUT:
% - OUTPUT directory is created
%         - creates a folder for icEEG and moves raw data sef file and resampled mat file into it
%         - creates a folder for hdEEG and moves raw data sef file and resampled mat file into it, as well as spatially resampled dbbEEG
% - icEEG: realigned intracranial EEG, channel labels and sampling frequency
% - hdEEG: realigned extracranial EEG, channel labels and sampling frequency

% OTHER
% - dual_clean_mkr:     make sure markers are cleaned, get rid of extra
%                       markers
%
% DEV
% - Maxime Baud, Feb 2017
% - refacto, Renaud Marquis, Aug 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%INPUTS
folder = fileparts(filename1);
fprintf('\nReading MicroMed icEEG data...\n\n')
[data1 des1] = dual_load_sef(filename1);
icchs = des1.channelnames';
fprintf('\nReading EGI hdEEG data...\n\n')
[data2 des2] = dual_load_sef(filename2);
hdchs = des2.channelnames';
load(mkrfile)
if ~isempty(mkr)
    
    start = mkr(1,1);
    stop = mkr(end,1);
    
    % Sampling frequ
    if nargin<4 || isempty(fs), fs=1000; end
    fs1 = median(diff(mkr(:,2)));
    fs2 = median(diff(mkr(:,3)));
    
    
    % FRAMES
    if nargin>6
        dc = mkr(:,1)<start | mkr(:,1)>stop;
        mkr(dc,:)=[];
    end
    nf = size(mkr,1); % number of frames
    % datapoints per frame
    dp1 = diff(mkr(:,2));
    dp2 = diff(mkr(:,3));
    
    
    % CLEAN CHANNELS
%     dc = cellfun(@isempty,icchs) | strcmp(icchs,'MKR1+ - '); % empty or marker+/- channels
    dc = cellfun(@isempty,icchs) | ~cellfun(@isempty,regexpi(icchs,'MKR')); % empty or marker+/- channels
    STIMmkr_labels = icchs(dc);
    STIMmkr_data = data1(dc,:);
    icchs(dc)=[];
    data1(dc,:)=[];
    dc = cellfun(@isempty,hdchs);
    hdchs(dc)=[];
    data2(dc,:)=[];
    
    % PRINT LOG
    cd(strrep(folder,'1_Raw','2_Resampled'))
    fid = fopen('Realignment_summary.txt','w'); % print log of what was realigned
    fprintf(fid,'Created: %s \nResampling frequency: %d\n\n',datestr(now,'dd-mmm-yyyy'),fs);
    fprintf(fid,'Frame \tMarkers \ticEEG [dp] \thdEEG [dp] \tInterpolation [dp]\n');
    
    % LOOP FRAMES
    icEEG = [];
    hdEEG = []; %preallocate
    for f=1:nf-1
        fprintf('%d ',f)
        if rem(f,20)==0 || f==nf-1, fprintf('\n'),end
        fl = round(fs*(mean([dp1(f)/fs1, dp2(f)/fs2]))); % frame length in dp
        X1 = mkr(f,2):mkr(f+1,2)-1;
        Y1 = data1(:,X1); % icEEG
        X2 = mkr(f,3):mkr(f+1,3)-1;
        Y2 = data2(:,X2); % hdEEG
        Xc = 1:fl;
        Y1b = interp1q(linspace(1,fl,dp1(f))',Y1',Xc')'; % interpolate icEEG
        if size(Y2,2)~=fl
            Y2b = interp1q(linspace(1,fl,dp2(f))',Y2',Xc')'; % interpolate hdEEG
        else
            Y2b = Y2;
        end
        
        % Concatenate EEGs
        icEEG = [icEEG Y1b];
        hdEEG = [hdEEG Y2b];
        
        fprintf(fid, '%d: \t%d to %d \t%d  \t%d \t%d\n',f, mkr(f,1), mkr(f+1,1),numel(X1),numel(X2),numel(Xc));  % print to log file
    end
    fclose(fid);
    
    if ~isempty(STIMmkr_data) && ~isempty(STIMmkr_labels)
        STIMmkr = [];
        for f=1:nf-1
        fprintf('%d ',f)
        if rem(f,20)==0 || f==nf-1, fprintf('\n'),end
        fl = round(fs*(mean([dp1(f)/fs1, dp2(f)/fs2]))); % frame length in dp
        X1 = mkr(f,2):mkr(f+1,2)-1;
        Y1 = STIMmkr_data(:,X1); % icEEG
        Xc = 1:fl;
        Y1b = interp1q(linspace(1,fl,dp1(f))',Y1',Xc')'; % interpolate icEEG
        
        % Concatenate MRK data
        STIMmkr = [STIMmkr Y1b];
        
        end
    end
    
    % STORAGE
    for n=1:2;
        switch n
            case 1; str = 'icEEG'; EEG = icEEG; labels = icchs; rawfile = filename1;
            case 2; str = 'hdEEG'; EEG = hdEEG; labels = hdchs; rawfile = filename2;
        end
        filepath = fullfile(strrep(folder,'1_Raw','2_Resampled'),str);
        if ~isdir(filepath),mkdir(strrep(folder,'1_Raw','2_Resampled'),str),end           % make folder if does not exist
        save(fullfile(filepath,str),'EEG', 'labels','fs')   % format that is accessible to Epitome
        
    end
    if ~isempty(STIMmkr_data) && ~isempty(STIMmkr_labels)
        save(fullfile(strrep(folder,'1_Raw','2_Resampled'),'STIM_MKR_chan.mat'),'STIMmkr_labels','STIMmkr','fs')
    end
    
else
    cd(strrep(folder,'1_Raw','2_Resampled'))
    fid = fopen('Realignment_summary.txt','w'); % print log of what was realigned
    fprintf(fid,'Created: %s \nERROR: EMPTY MARKERS!\n\n',datestr(now,'dd-mmm-yyyy'));
    fclose(fid);
end

end


