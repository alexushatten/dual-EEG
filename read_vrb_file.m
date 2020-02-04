function VerboseFile = read_vrb_file(Filename)
% read_vrb_file: read Cartool verbose (.vrb) file, useful for finding info
% in processed files using regex
%
% VerboseFile = read_vrb_file(Filename)
%
%  Inputs
% --------
% Filename: path to .vrb file 
%
%  Outputs
% ---------
% VerboseFile: cell array of strings with the content of the text file
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, 2018
%-------------------------------------------------------------------------

formatSpec = '%s%[^\n\r]';
fileID = fopen(Filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '',  'ReturnOnError', false);
dataArray{1} = strtrim(dataArray{1});
fclose(fileID);
VerboseFile = dataArray{:, 1};

end
