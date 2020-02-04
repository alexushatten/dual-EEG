function [data, numchannels, numtimeframes, samplingrate]=read_eph(ephfilepath)
% openeph: opens a Cartool evoked potential data file (.ep(h))
%
% Usage:
%
% [data, numchannels, numtimeframes, samplingrate]=read_eph(ephfilepath)
%
% Inputs: full path and name of the file to open
%
% Outputs: data as a 2-D numeric array where dimension 1 contains the
% timeframes, dimension 2 contains the channels; number of channels, number
% of timeframes (and sampling rate for .eph files) as 1-D numeric arrays
%
% Cartool: http://brainmapping.unige.ch/Cartool.htm
%
% author of this script: pierre.megevand@medecine.unige.ch
% refacto, Renaud Marquis @FBMlab, June 2018


% open filename for reading in text mode
fid=fopen(ephfilepath,'rt');

% for .eph files
if strcmp(ephfilepath(end-3:end),'.eph')==1
    
    % read header
    eph=textscan(fid,'%s','delimiter','/n');
    eph=eph{1};
    numchannels=sscanf(eph{1},'%f',1);
    numtimeframes=sscanf(eph{1},'%*f %f',1);
    samplingrate=sscanf(eph{1},'%*f %*f %f',1);
    
    % prepare for reading data
    formatstring='%f';
    if numchannels>1
        for i=1:numchannels-1
            formatstring=[formatstring ' %f'];
        end
    end

    % read data
    if size(eph,1)==numtimeframes+1
        data=zeros(numtimeframes,numchannels);
        for i=1:numtimeframes
            data(i,:)=sscanf(eph{i+1},formatstring);
        end
    elseif size(eph,1)==numtimeframes % this can happen e.g. with EEG epochs concatenated after averaging
        data=zeros(numtimeframes-1,numchannels);
        for i=1:numtimeframes-1
            data(i,:)=sscanf(eph{i+1},formatstring);
        end
    end

% for .ep files
elseif strcmp(ephfilepath(end-2:end),'.ep')==1
    data='';
    while ~feof(fid)
        thedataline=fgetl(fid);
        data=strvcat(data,thedataline);
    end
    data=str2num(data);
    numtimeframes=size(data,1);
    numchannels=size(data,2);
    samplingrate=0;
    
else
    error('incorrect file type');
end

% close file
fclose(fid);

end