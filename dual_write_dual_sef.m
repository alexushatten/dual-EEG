function dual_write_dual_sef( savefilename, icEEG, hdEEG, samplingrate, WriteGFP, icEEGchan, hdEEGchan )
% dual_write_dual_sef: writes .sef file with simultaneous scalp and
% intracranial EEG data that will be scaled such that they can be
% conveniently viewed simultaneously in Cartool
%
% dual_write_dual_sef( savefilename, icEEG, hdEEG, samplingrate, icEEGchan, hdEEGchan )
%
%  Inputs
% --------
% savefilename: full path and name of the file to save
% icEEG: timeframes x channels array of intracranial EEG
% hdEEG: timeframes x channels array of scalp EEG
% samplingrate: integer, sampling rate of icEEG and hdEEG (should be equal!)
% WriteGFP: logical, whether to write GFP as an additional channel
%           [default = true]
% icEEGchan (optional): channel names of icEEG as cell array of strings
% hdEEGchan (optional): channel names of hdEEG as cell array of strings
%
%  Outputs
% ---------
% separate .sef files with scaled icEEG and hdEEG traces, but scaled such
% that they can be easily visualized simultaneously in Cartool
%
%-------------------------------------------------------------------------
% Nota bene: bad channels are assumed to be removed, if it is not the case,
% visualization will be bad
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, September 2018
%-------------------------------------------------------------------------

if nargin < 5
    WriteGFP = true;
end

if ~sum(size(icEEG)==size(hdEEG))
    error('icEEG & hdEEG size mismatch')
end

TimeframesDim = find(size(icEEG)==size(hdEEG));

if isempty(TimeframesDim)
    error('Input dimension does not fulfill requirements, matching dimension not found, aborting...')
elseif TimeframesDim ~= 1
    warning('Input dimension does not fulfill requirements but matching dimension found, transposing input traces')
    icEEG = icEEG';
    hdEEG = hdEEG';
end

NchanIC = size(icEEG,2);
NchanHD = size(hdEEG,2);
TotalNchan = NchanIC+NchanHD;

% scale data but preserve polarity
scaledICeeg = icEEG/range(icEEG(:));
scaledHDeeg = hdEEG/range(hdEEG(:));

% make fake electrodes to match number of electrodes
if ~WriteGFP
    scaledICeeg = [scaledICeeg,zeros(size(scaledICeeg,1),TotalNchan-NchanIC)];
    scaledHDeeg = [scaledHDeeg,zeros(size(scaledHDeeg,1),TotalNchan-NchanHD)];
    if nargin > 5
        hdEEGchan(end+1:TotalNchan) = {'null'};
        icEEGchan(end+1:TotalNchan) = {'null'};
    end
elseif WriteGFP
    GFP = computegfp(scaledHDeeg);
    scaledICeeg = [scaledICeeg,GFP,zeros(size(scaledICeeg,1),TotalNchan-NchanIC)];
    scaledHDeeg = [scaledHDeeg,zeros(size(scaledHDeeg,1),TotalNchan-NchanHD+1)];
    if nargin > 5
%         hdEEGchan(end+1) = {'GFP'};
        hdEEGchan(end+1:TotalNchan+1) = {'null'};
        icEEGchan(end+1) = {'GFP'};
%         icEEGchan(end+1:TotalNchan+1) = {'null'};
        icEEGchan(end+1:TotalNchan+1) = {'null'};
    end
else 
    error('Check WriteGFP input arg format!')
end

if nargin > 5
    write_sef(spm_file(savefilename,'suffix','_icEEG'),scaledICeeg,samplingrate,icEEGchan);
    write_sef(spm_file(savefilename,'suffix','_hdEEG'),scaledHDeeg,samplingrate,hdEEGchan);
elseif nargin < 6
    write_sef(spm_file(savefilename,'suffix','_icEEG'),scaledICeeg,samplingrate);
    write_sef(spm_file(savefilename,'suffix','_hdEEG'),scaledHDeeg,samplingrate);
else
    error('Check number of input arguments')
end

end
