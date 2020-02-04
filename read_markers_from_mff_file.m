function [RelativeTime, Code, Label, Description, Duration, RecordStartTime] = read_markers_from_mff_file(FilePath)
% read_markers_from_mff_file  reads markers from EGI .mff datafile
%
% USAGE
% [RelativeTime, Code, Label, Description, Duration, RecordStartTime] = read_markers_from_mff_file(FilePath)
%
% INPUTS
% FilePath: string, path to .mff file
%
% OUTPUTS
% RelativeTime: double, time stamps of markers in milliseconds
% Code: cell of strings, marker codes
% Label: cell of strings, marker labels
% Description: cell of strings, marker descriptions (usually not useful)
% Duration: double, marker duration (in microseconds, usually not useful)
% RecordStartTime: string, time stamp indicating starting time of recording
% (usually not useful)
%
% Nota bene: looks for time markers automatically generated from experimental
% setup (triggers, or markers that cannot be edited in NetStation) in
% dualEEG experiments matching the keyword 'DIN' (which stands possibly for Digital INput)
%
% R Marquis, 2017-07-11

warning('Depending on EGI''s NetStation version, XML structure might change slightly and output variables might be unreliable... So always check the output! (and update children''s positions (around lines 88-112) if not doable...)')

if nargout > 5
    warning('Recording starting time in RecordStartTime might be unreliable')
end

PkgFiles = dir(strcat(FilePath,filesep,'*.xml'));

% if info.xml file not found in .mff file, start time of recording cannot
% be inferred:
if exist(strcat(FilePath,filesep,'info.xml'))~=2
    error('"info.xml" not found in %s, cannot proceed further... :-(',FilePath)
else
    InfoXML = strcat(FilePath,filesep,'info.xml');
end

% Events seem to be always stored in an .xml file starting with
% "Events_", unless someone manually changed the content of the .mff file.

IdxPossibleEventsFile = find(~cellfun(@isempty,regexp({PkgFiles.name},'^Events_*')));

if isempty(IdxPossibleEventsFile)
    error('"Events_***.xml" not found in %s, cannot proceed further... :-(',FilePath)
end

% Any filename can be found based on the name defined by the user when
% creating the track (e.g. 'Events_User Markup.xml', 'Events_newTrack.xml',
% 'Events_selection.xml').
% However here we know that the events we are looking for
% should match the "DIN" keyword because they are triggers that were
% generated automatically by the experimental setup during recording.

% Most of the times, it is named 'Events_255 DINs.xml' but in one case it
% is named 'Events_DIN_1.xml'.

EventsFileWithinPossible = find(~cellfun(@isempty,regexp({PkgFiles(IdxPossibleEventsFile).name},'DIN')));

if isempty(EventsFileWithinPossible)
    error('"Events_***.xml" file found in %s, but digital inputs are missing (no match for keyword "DIN"): try another keyword for User Markup events, etc...',FilePath)
end

EventsXML = strcat(FilePath,filesep,PkgFiles(IdxPossibleEventsFile(EventsFileWithinPossible)).name);


Events = parseXML(EventsXML);

% seems that childrens with name "event" are more relevant than others...
IdxEvent = zeros(size(Events.Children),'single')';
for ntemp = 1:size(Events.Children,2)
    IdxEvent(ntemp) = strcmp(Events.Children(ntemp).Name,'event');
end
IdxEvents = find(IdxEvent);

% seems that childrens of childrens at positions 2, 4, 6, 8 & 10 are more
% relevant than others...
BeginTime = cellstr(repmat('',size(IdxEvents)));
Duration = zeros(size(IdxEvents));
Code = cell(size(IdxEvents));
Label = cell(size(IdxEvents));
Description = cell(size(IdxEvents));
for n = 1:length(IdxEvents)
    % lazy coding here: I assume here that the order of the different fields
    % remains the same to avoid looping to match the name of the fields...
    if isfield(Events.Children(IdxEvents(n)).Children(2).Children,'Data')
        BeginTime{n} = Events.Children(IdxEvents(n)).Children(2).Children.Data;
    else
        error('Could not find markers onset, cannot compute relative time')
    end
    if isfield(Events.Children(IdxEvents(n)).Children(4).Children,'Data')
        Duration(n) = str2double(Events.Children(IdxEvents(n)).Children(4).Children.Data); % in 탎econds
    else
        Duration(n) = nan;
    end
    if isfield(Events.Children(IdxEvents(n)).Children(6).Children,'Data')
        Code{n} = Events.Children(IdxEvents(n)).Children(6).Children.Data;
    else
        Code{n} = '';
    end
    if isfield(Events.Children(IdxEvents(n)).Children(8).Children,'Data')
        Label{n} = Events.Children(IdxEvents(n)).Children(8).Children.Data;
    else
        Label{n} = '';
    end
    if isfield(Events.Children(IdxEvents(n)).Children(10).Children,'Data')
        Description{n} = Events.Children(IdxEvents(n)).Children(10).Children.Data;
    else
        Description{n} = '';
    end    
end

% time stamps of markers are in absolute date, so one needs:
%   a) the starting date of the recording, found in 'info.xml':
StartTime = parseXML(InfoXML);
StartTimeTag = zeros(size(StartTime.Children),'single')';
for ntemp = 1:size(StartTime.Children,2)
    StartTimeTag(ntemp) = strcmp(StartTime.Children(ntemp).Name,'recordTime');
end

RecordStartTime = StartTime.Children(find(StartTimeTag)).Children.Data;

%   b) to convert absolute time stamps to relative time stamps
RelativeTime = zeros(size(IdxEvents));
for n = 1:length(IdxEvents)
    [~,~,RelativeTime(n)] = get_time_vector_from_mff_file(BeginTime{n},RecordStartTime);
end

end

%------------- INTERNAL FUNCTIONS -------------

function [TimeVector, StartingTime, TimeDiff] = get_time_vector_from_mff_file(TimeStamp,StartTime)
% TimeVector  Get time vector (datevec) from EGI .mff datafile timestamp
% 
% INPUTS
% TimeStamp: string, time stamp with the format yyyy-mm-ddThh:MM:SS.도도도+/-GMT
% StartTime [optional]: string, starting time of recording, same format as TimeStamp
%
% OUTPUTS
% TimeVector: double, time stamp with:
%    1st element: years 
%    2nd element: months
%    3rd element: days
%    4th element: hours
%    5th element: minutes
%    6th element: seconds
%    7th element: milliseconds
%
% StartingTime: double, starting time with format as TimeVector
% TimeDiff: double, time difference between TimeStamp and StartTime
% 
% Nota bene: GMT+/-xxx is simply ignored
%
% R Marquis, 2017-07-11

if nargin < 2
    StartingTime = 0;
else
    StartingTime = datevec(datenum(StartTime(1:19),'yyyy-mm-ddTHH:MM:SS'));
    StartingTime(end+1) = str2double(StartTime(21:26)); % for microseconds, because Matlab only handles milliseconds
end

% datetime(strcat(datestr(datenum(RawRecordStartTime(1:19),'yyyy-mm-ddTHH:MM:SS')),'.',RawRecordStartTime(21:26)),'Format','yyyy-MM-dd,hh:mm:ss.SSSSSS')
TimeVector = datevec(datenum(TimeStamp(1:19),'yyyy-mm-ddTHH:MM:SS'));
TimeVector(end+1) = str2double(TimeStamp(21:26)); % for microseconds, because Matlab only handles milliseconds

% TimeDiff = TimeVector-StartingTime;
TimeDiffTMP = etime(TimeVector(1:6),StartingTime(1:6)); % in seconds

% deal with 탎

if StartingTime(7)>TimeVector(7)
    TimeDiffmicroSec = 1e6+(TimeVector(7)-StartingTime(7));
    TimeDiffTMP = TimeDiffTMP-1;
else
    TimeDiffmicroSec = (TimeVector(7)-StartingTime(7));
end

TimeDiff = round((TimeDiffTMP*1000000+TimeDiffmicroSec)/1000); % convert to 탎 and then back to ms

end
