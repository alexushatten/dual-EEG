function [icEEG,icchs,hdEEG,hdchs,fs] = dual_realign(icEEGfilename, hdEEGfilename, merged_marker_filename, resampf)  %#RM, 16.06.17 09:00:19: insert [filename1, filename2, fs, mkrfile ,start, stop] as input arguments

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
% - filename1, data 1:     intracranial data (micromed 2048Hz)
% - filename2, data 2:     scalp data (EGI 1000 Hz)
% - marker file:    markers with labels (col 1), and datapoints for icEEG (col 2) and hdEEG (col 3)
% - start/stop: optional, default 1:end, markers label at which to start
%               and stop
% - resampf:    optional, resampling frequency, 1000 Hz by default
%
% OUTPUT: 
% - icEEG: realigned intracranial EEG 
% - hdEEG: realigned extracranial EEG
% - chs: channel labels ic-... for intracranial and hd-... for extracranial


% OTHER
% - dual_clean_mkr:     make sure markers are cleaned, get rid of extra
%                       markers
%
% DEV
% - Maxime Baud, Feb 2017
% - Renaud Marquis, refacto, July 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% OUTPUTS
EEG = [];

% PROMPT inputs
filename1 = icEEGfilename;
% filename1 = input('ic-EEG file: ','s'); %#RM, 04.07.17 08:44:34
[data1 des1] = dual_load_sef(filename1);
chs1 = des1.channelnames';
filename2 = hdEEGfilename;
% filename2 = input('HD-EEG file: ','s'); %#RM, 04.07.17 08:44:42
[data2 des2] = dual_load_sef(filename2);
chs2 = des2.channelnames';
fs = resampf;
% fs = input('Resampling frequency: 1000'); %#RM, 04.07.17 08:45:58
if isempty(fs), fs=1000; end
mkrfile = merged_marker_filename;
% mkrfile = input('Merged markers file: ','s'); %#RM, 04.07.17 08:48:33
load(mkrfile)
start = input(sprintf('Start at marker: %d', mkr(1,1)),'s');
stop = input(sprintf('Stop at marker: %d', mkr(end,1)),'s');

%#RM, 04.07.17 18:10:48
folder = fileparts(filename1);
folder2 = fileparts(filename2);
if ~strcmp(folder,folder2)
    warning('hdEEG and icEEG not in the same folder, log file will be saved at icEEG file folder, i.e. %s',folder)
end
if isempty(folder)
    folder = cd;
end
%#RM, 04.07.17 18:10:48

% Sampling frequ
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


% CLEAN icCHANNELS
dc = cellfun(@isempty,chs1) | strcmp(chs1,'MKR1+ - '); % empty or marker+/- channels
chs1(dc)=[];
data1(dc,:)=[];
dc = cellfun(@isempty,chs2);
chs2(dc)=[];
data2(dc,:)=[];
chs = [chs1;chs2];

% PRINT LOG
cd(folder)
fid = fopen('Realignment_summary.txt','w'); % print log of what was realigned
fprintf(fid,'Created: %s \nResampling frequency: %d\n\n',datestr(now,'dd-mmm-yyyy'),fs);
fprintf(fid,'Frame \tMarkers \ticEEG [dp] \thdEEG [dp] \tInterpolation [dp]\n');

% LOOP FRAMES
for f=1:nf-1
    fprintf('%d ',f)
    if rem(f,20)==0 || f==nf-1, fprintf('\n'),end
    fl = round(fs*(mean([dp1(f)/fs1, dp2(f)/fs2]))); % frame length in dp
    X1 = mkr(f,2):mkr(f+1,2)-1;
    Y1 = data1(:,X1); % icEEG
    X2 = mkr(f,3):mkr(f+1,3)-1;
    Y2 = data2(:,X2); % hdEEG
    Xc = 1:fl;
    Y1b = interp1q(linspace(1,fl,dp1(f))',Y1',Xc')'; % interpolate icEEG %#RM, 16.06.17 09:04:27: why not interp1 instead? (should be faster AFAIK)
    if size(Y2,2)~=fl 
        Y2b = interp1q(linspace(1,fl,dp2(f))',Y2',Xc')'; % interpolate hdEEG
    else
        Y2b = Y2;
    end
   
    Yc = [Y1b;Y2b];
    EEG = [EEG Yc]; % concatenate EEG
    fprintf(fid, '%d: \t%d to %d \t%d  \t%d \t%d\n',f, mkr(f,1), mkr(f+1,1),numel(X1),numel(X2),numel(Xc));  % print to log file
end
fclose(fid);
if nargout>2
    ic_idx = numel(chs1);
    icEEG = EEG(1:ic_idx,:);
    hdEEG = EEG(ic_idx+1:end,:);
    icchs = chs(1:ic_idx);
    hdchs = chs(ic_idx+1:end);
end