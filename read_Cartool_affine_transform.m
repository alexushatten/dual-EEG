function A = read_Cartool_affine_transform(filename)
%read_Cartool_affine_transform: read .txt output file from Cartool after
%inverse solution computation
%
% Usage:
%-------------------------------------------------------------------------
% 
% A = read_Cartool_affine_transform(filename)
%
% Inputs:
%-------------------------------------------------------------------------
%
% filename: string, full file path
%
% Outputs:
%-------------------------------------------------------------------------
%
% A: affine transformation matrix
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, March 2018
%-------------------------------------------------------------------------

delimiter = '\t';
if nargin<=2
    startRow = 1;
    endRow = inf;
end

formatSpec = '%f%f%f%f%[^\n\r]';

fileID = fopen(filename,'r');

dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

fclose(fileID);

A = [dataArray{1:end-1}];

end

