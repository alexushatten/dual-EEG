function dual_concatenate_sef( FilesArray, OutputDir, Check )
% dual_concatenate_sef: concatenate .sef files (handles .mrk files as well)
%
% dual_concatenate_sef( FilesArray, OutputDir, Check )
%
%  Inputs
% --------
% FilesArray: cell array of strings with file paths to .sef files
% OutputDir: string, directory for output (default = directory of first
% file)
% Check: whether to check sampling frequency and number of
% channels across files to concatenate (default = true)
%
%  Outputs
% ---------
% SEF file where data was concatenated, with associated MRK files (filename
% will be concatenated input filenames).
%
%-------------------------------------------------------------------------
% NB: as long as each individual .sef file holds in the computer's memory,
% the resulting file with concatenated EEG data will not make the
% computer's memory explode because it will read and write parts after
% parts ;-)
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, October 2018
%-------------------------------------------------------------------------

if nargin < 3
    Check = true;
elseif nargin < 2
    OutputDir = fileparts(FilesArray{1});
end

% append .sef
% DataAll = [];

Hdr = dual_load_sef_hdr(FilesArray{1});
% open savefilename for writing
fid=fopen(fullfile(OutputDir,[strtrim(vectorize(char(spm_file(FilesArray(:),'basename'))')'),'.sef']),'w');

% get total number of samples across all files
NumSamples = zeros(length(FilesArray),1);
for f = 1:length(FilesArray)
    Hdr = dual_load_sef_hdr(FilesArray{f});
    if Check && f>1
        if Hdr.samplingfreq~=HdrPrev.samplingfreq
            fclose(fid);
            error('Sampling frequency is not the same in %s and %s',FilesArray{f},FilesArray{f-1})
        elseif Hdr.nbchannels~=HdrPrev.nbchannels
            fclose(fid);
            error('Number of channels is not the same in %s and %s',FilesArray{f},FilesArray{f-1})
        end
    end
    HdrPrev = Hdr;
    NumSamples(f) = Hdr.nbsamples;
end

%write fixed part of header
fwrite(fid,'SE01','int8');
fwrite(fid,Hdr.nbchannels,'int32');
fwrite(fid,Hdr.Naux,'int32');
fwrite(fid,sum(NumSamples),'int32');
fwrite(fid,Hdr.samplingfreq,'float32');
fwrite(fid,Hdr.year,'int16');
fwrite(fid,Hdr.month,'int16');
fwrite(fid,Hdr.day,'int16');
fwrite(fid,Hdr.hour,'int16');
fwrite(fid,Hdr.minute,'int16');
fwrite(fid,Hdr.second,'int16');
fwrite(fid,Hdr.msecond,'int16');

% define and write variable part of header
for i=1:Hdr.nbchannels
    currentchannel=uint8(Hdr.channelnames{i});
    for j=1:size(currentchannel,2)
        fwrite(fid,currentchannel(j),'int8');
    end
    for k=j+1:8
        fwrite(fid,0,'int8');
    end
end

% write the data:
for f = 1:length(FilesArray)
    [Data,~] = dual_load_sef(FilesArray{f});
    fwrite(fid,Data,'float32');
end

% close file
fclose(fid);

%% check for presence of .sef.mrk files, if such, concatenate them as well
MrkThere = false(length(FilesArray),1);
for f = 1:length(FilesArray)
    MrkPath{f} = spm_file(FilesArray{f},'ext','sef.mrk'); %#ok<AGROW>
    MrkThere(f) = exist(MrkPath{f},'file')==2;
end
MrkPath = MrkPath';
MrkPath = MrkPath(MrkThere);
CumNumSamples = cumsum(NumSamples);
FileMrgMrk = CumNumSamples+1;
FileMrgMrk = FileMrgMrk(1:end-1);
StartTime = 0;
MarkerStart = [];
MarkerEnd = [];
MarkerLabel = [];
for f = 1:length(MrkPath)
    try
        [tS, tE, tL] = read_mrk_file(MrkPath{f});
    catch ME
        warning([ME.message,' Trying fallback marker reading function instead.'])
        [tS, tE, tL] = read_mrk_Cartool(MrkPath{f});
    end
    MarkerStart = [MarkerStart; StartTime+tS]; %#ok<AGROW>
    MarkerEnd = [MarkerEnd; StartTime+tE]; %#ok<AGROW>
    MarkerLabel = [MarkerLabel; tL]; %#ok<AGROW>
    StartTime = StartTime+NumSamples(f);
end

% add markers where files were merged
FinalMarkerStart = [MarkerStart;FileMrgMrk];
FinalMarkerEnd = [MarkerEnd;FileMrgMrk];
FinalMarkerLabel = [MarkerLabel;cellstr(repmat('"FileMergePoint"',length(FileMrgMrk),1))];
[~,IdxSort] = sort(FinalMarkerStart);

if max(FinalMarkerStart)>CumNumSamples(end) || max(FinalMarkerEnd)>CumNumSamples(end)
    warning('Some markers are out-of-range!')
end

if ~isempty(MrkPath)
    write_mrk_file_Cartool(fullfile(OutputDir,[strtrim(vectorize(char(spm_file(FilesArray(:),'basename'))')'),'.sef.mrk']),FinalMarkerStart(IdxSort),FinalMarkerEnd(IdxSort),FinalMarkerLabel(IdxSort));
end
% write_sef(fullfile(OutputDir,[strtrim(vectorize(char(spm_file(FilesArray(:),'basename'))')'),'.sef']),DataAll',Hdr.samplingfreq,Hdr.channelnames);

end

