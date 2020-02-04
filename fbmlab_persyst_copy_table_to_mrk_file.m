function fbmlab_persyst_copy_table_to_mrk_file(Filename, SamplingRate)
% fbmlab_persyst_copy_table_to_mrk_file: writes .mrk file based on csv file
% containing PERSYST Reveal results (copied as table to clipboard)
%
% fbmlab_persyst_copy_table_to_mrk_file(Filename,SamplingRate)
%
%  Inputs
% --------
% Filename: string, filepath to .csv file (comma (or semicolon) -separated
% values
% SamplingRate : integer, number of samples per second for the requested
% .mrk
%
%  Outputs
% ---------
% .mrk file with first column of .csv file as event (marker) label, onset and
% offset of events.
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, January 2019
%-------------------------------------------------------------------------

%% Initialize variables.
delimiter = ',';
% if nargin<=2
startRow = 2;
endRow = inf;
% end

%% Format string for each line of text:
%   column1: text (%s)
%	column2: text (%s)
%   column3: double (%f)
%	column4: text (%s)
%   column5: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%f%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(Filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
Text = dataArray{:, 1};
Time = dataArray{:, 2};

if all(cellfun(@isempty,Time))
    warning('Microsoft, please, .csv means COMMA-separated values, not semicolon separated values...')
    warning('Okay, I am gonna deal with this anyway.')

    delimiter = ';';

    fileID = fopen(Filename,'r');
    dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
    for block=2:length(startRow)
        frewind(fileID);
        dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
        for col=1:length(dataArray)
            dataArray{col} = [dataArray{col};dataArrayBlock{col}];
        end
    end
    fclose(fileID);
    
    Text = dataArray{:, 1};
    Time = dataArray{:, 2};
    
end

TimeChar = char(Time);
Time = cellstr(TimeChar(:,4:end));
TempTime = regexp(Time,':|\.','split');

Hours = nan(length(TempTime),1);
Minutes = nan(length(TempTime),1);
Seconds = nan(length(TempTime),1);
MilliSeconds = nan(length(TempTime),1);
for tt = 1:length(TempTime)
    
    Hours(tt) = str2num(TempTime{tt}{1});
    Minutes(tt) = str2num(TempTime{tt}{2});
    Seconds(tt) = str2num(TempTime{tt}{3});
    MilliSeconds(tt) = str2num(TempTime{tt}{4});
end

Duration = dataArray{:, 3};
% Origin = dataArray{:, 4};
% Sort = dataArray{:, 5};

T1 = datenum(Hours,Minutes,Seconds)*SamplingRate + round(MilliSeconds/1000*SamplingRate);
T2 = Duration*SamplingRate + T1;

write_mrk_file_Cartool(spm_file(Filename,'ext','mrk'),T1,T2,Text)

end
