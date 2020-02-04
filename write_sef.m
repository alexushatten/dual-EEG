function write_sef(savefilename,thedata,samplingrate,varargin)
% write_sef: saves data as a Cartool simple EEG data file (.sef)
%
% write_sef(savefilename,thedata,samplingrate,channelnames,header)
%
%  Inputs
% --------
% savefilename: full path and name of the file to save
% thedata: [time x channels] EEG traces as a 2-D numeric array, where
% dimension 1 contains the timeframes, dimension 2 contains the channels
% samplingrate: samplingrate as 1-D numeric array
% channelnames (optional): channel names as cell array of strings
% header (optional): header for absolute time information, following
%                    this format:
%                               header.year:        int16, year (of acquisition)
%                               header.month:       int16, month
%                               header.day:         int16, day
%                               header.hour:        int16, hour
%                               header.minute:      int16, minute
%                               header.second:      int16, second
%                               header.msecond:     int16, millisecond
%
%  Outputs
% ---------
% .sef file
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Cartool: http://brainmapping.unige.ch/Cartool.htm
%
% Based on savesef.m by pierre.megevand@medecine.unige.ch
% added channel names, renaud.marquis@unige.ch, 2017-08-22
% fixed empty channel names, Renaud Marquis @ FBM lab, December 2018
% #RM@FBMlab: added optional header information (absolute time)
%
% Renaud Marquis @ FBMlab, June 2019
%-------------------------------------------------------------------------

% define fixed part of header
version='SE01';
numchannels=size(thedata,2);
numauxchannels=0;
numtimeframes=size(thedata,1);
if nargin > 4
    year = varargin{2}.year;
    month = varargin{2}.month;
    day = varargin{2}.day;
    hour = varargin{2}.hour;
    minute = varargin{2}.minute;
    second = varargin{2}.second;
    millisecond = varargin{2}.msecond;
else
    year=0;
    month=0;
    day=0;
    hour=0;
    minute=0;
    second=0;
    millisecond=0;
end

% open savefilename for writing
fid=fopen(savefilename,'w');

%write fixed part of header
fwrite(fid,version,'int8');
fwrite(fid,numchannels,'int32');
fwrite(fid,numauxchannels,'int32');
fwrite(fid,numtimeframes,'int32');
fwrite(fid,samplingrate,'float32');
fwrite(fid,year,'int16');
fwrite(fid,month,'int16');
fwrite(fid,day,'int16');
fwrite(fid,hour,'int16');
fwrite(fid,minute,'int16');
fwrite(fid,second,'int16');
fwrite(fid,millisecond,'int16');

% define and write variable part of header
if nargin<4
    for i=1:numchannels
        fwrite(fid,101,'int8');
        currentchannel=uint8(num2str(i));
        for j=1:size(currentchannel,2)
            fwrite(fid,currentchannel(j),'int8');
        end
        for k=j+2:8
            fwrite(fid,0,'int8');
        end
    end
else
    if numel(varargin{1})~=size(thedata,2)
        fclose all;
        error('Mismatch between input channel names and number of channels in EEG traces')
    end
    for i=1:numchannels
%         TempCurChan = regexprep(varargin{1}{i},'-|+','');
        TempCurChan = varargin{1}{i};
        if isempty(TempCurChan)
%             TempCurChan = 'none'; % because Cartool does not like channels with empty labels...
            TempCurChan = ' '; % because Cartool does not like channels with empty labels...
        end
        currentchannel=uint8(TempCurChan);
        for j=1:size(currentchannel,2)
            fwrite(fid,currentchannel(j),'int8');
        end
        for k=j+1:8
            fwrite(fid,0,'int8');
        end
    end
end

% write data
fwrite(fid,thedata','float32');

% close file
fclose(fid);

end