% This scirpt converts csv file in relative time from Persyst in mrk file in time frames for Cartool
% Vincent ROCHAS - FBMLab 2018.11 - matlab R2016a
% Inputs:   the csv file path and file name.
%           the sampling frequency.
%
% #RM@FBMlab: fixed issue with input args, 2019-01-29

function fbmlab_csv_to_mrk_file(varargin)
%% Initialize variables.
if nargin>1 && ~isempty(varargin{1})
    filepath=varargin{1};
    filename=varargin{2};
    file = fullfile(filepath,filename); % #RM@FBMlab: fixed issue with input args
else
    [filename,filepath]=uigetfile('*.csv','Select a .csv file');
    file=[filepath,filename];
end

if nargin>2 && ~isempty(varargin{1})
    fs=varargin{1};
else
    dlg_title = 'Enter a sample frequency';
    prompt = {'Enter a sample frequency'};
    num_lines = 1;  def = {'256'};
    options.Resize='on';
    fs = inputdlg(prompt,dlg_title,num_lines,def,options);
    fs = str2double(fs);
end

delimiter = ',';
startRow = 2;

%% Format string for each line of text:
%   column1: text (%q)
%	column2: text (%q)
%   column3: double (%f)
%	column4: text (%q)
%   column5: text (%q)
%	column6: text (%q)
%   column7: text (%q)
%	column8: text (%q)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%q%q%f%q%q%q%q%q%[^\n\r]';

%% Open the text file.
fileID = fopen(file,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
Time1 = dataArray{:, 1};
Channel1 = dataArray{:, 2};
Perception1 = dataArray{:, 3};
Sign1 = dataArray{:, 4};
Duration1 = dataArray{:, 5};
Height1 = dataArray{:, 6};
Angle1 = dataArray{:, 7};
Group1 = dataArray{:, 8};


%% Clear temporary variables
clearvars delimiter startRow formatSpec fileID dataArray ans;

%% Collect time in time frame and collet trigger names in a table
table=cell(size(Time1,1),3);	% A cell array with 3 columns for Cartool -> 1:TF 2:TF 3:marker name
for m = 1:size(Time1,1) 
% TimeString = '00:00:05.720'; % the time string from the csv file
TimeString = char(Time1(m)); % the time string from the csv file
TimeString = TimeString(4:end); % remove the day info
t1 = datevec(TimeString,'HH:MM:SS'); % separate into hours, minutes and seconds
ZeroString = '00:00:00';    % a string of t0
t0 = datevec(ZeroString,'HH:MM:SS');	% separate into hours, minutes and seconds
t1seconds = etime(t1,t0);	% the diff in seconds of t1 from 0
% t1seconds=t1seconds+str2num(TimeString(end-3:end)); % plus the millisecond part that was missing
t1seconds=t1seconds+str2double(TimeString(end-3:end)); % plus the millisecond part that was missing
t1tf=t1seconds*fs;  % convert second to time frame
table{m,1}=t1tf;
table{m,2}=t1tf;

mrkname=char(Channel1(m)); % the name string from the csv file
table{m,3}=mrkname;
end

% Create the marker file
fileout= [file(1:end-4) '.mrk'];
fid_mrk=fopen(fileout,'wt');
% Fill in the empty file
fprintf(fid_mrk, '%s \n', 'TL02');    % first line with TL02 for Cartool
for k=1:size(Time1,1) 
    fprintf(fid_mrk, '%d \t %d \t %s \n', round(table{k,1}), round(table{k,2}), ['"' table{k,3} '"']);
end
% Close it
fclose(fid_mrk);
clear all
end
