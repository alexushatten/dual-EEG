function Montage = read_mtg_file(Filename)
% read_mtg_file: read Cartool montage (.mtg) file
%
% Montage = read_mtg_file(Filename)
%
%  Inputs
% --------
% Filename: path to .vrb file 
%
%  Outputs
% ---------
% Montage: [n x 2] cell array of strings with pairs of channel names
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, February 2019
%-------------------------------------------------------------------------

formatSpec = '%s%[^\n\r]';
fileID = fopen(Filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '\t', 'WhiteSpace', '',  'ReturnOnError', false);
dataArray{1} = strtrim(dataArray{1});
fclose(fileID);
if size(dataArray,1)~=1 || size(dataArray,2)~=2
    error('Improperly formatted .mtg file')
else
    Montage = dataArray{:,1};
    Montage(:,2) = dataArray{:,2};
    Montage = Montage(2:end,:);
end
end
