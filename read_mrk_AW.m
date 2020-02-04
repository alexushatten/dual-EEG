function [Label,Code,Onset,Duration,Comment,Color] = read_mrk_AW(Filename)
% read_mrk_AW: read AnyWave .mrk file
%
% [MarkerLabel,Code,Onset,Duration] = read_mrk_AW(Filename)
%
%  Inputs
% --------
% Filename: path to AnyWave .mrk file
%
%  Outputs
% ---------
% The variable name says all!
% Please not that Onset and Duration are in seconds!
%
% http://meg.univ-amu.fr/wiki/AnyWave:ADES#The_marker_file_.28.mrk.29
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, March 2019
%-------------------------------------------------------------------------

%% Initialize variables.
delimiter = '\t';
startRow = 2;
endRow = inf;

%% Format string for each line of text:
%   column1: text (%q)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: text (%q)
%	column6: text (%q)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%q%f%f%f%q%q%[^\n\r]';

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
Label = dataArray{:, 1};
Code = dataArray{:, 2};
Onset = dataArray{:, 3};
Duration = dataArray{:, 4};
Comment = dataArray{:, 5};
Color = dataArray{:, 6};

if ~(sum(cellfun(@isempty,regexp(Comment,'^#'))) ==  length(Label))
    warning('Parsing of .mrk file might have possibly go wrong, so unless some targets (comments) start with "#", the parsing mixed colors and targets')
end
if sum(cellfun(@isempty,regexp(Color,'^#')) - cellfun(@isempty,Color))>0
    warning('Parsing of .mrk file might have possibly go wrong, because some color codes do not start with "#"')
end

end
