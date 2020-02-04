function [x,y,z,name,Xyz] = read_xyz_file(Filename)
%read_xyz_file: reads Cartool .xyz file
%
% Usage:
%-------------------------------------------------------------------------
% 
% [x,y,z,name,Xyz] = read_xyz_file(Filename)
%
% Outputs:
%-------------------------------------------------------------------------
%
% x, y, z coordinates
%
% name: name of electrode / contact
%
% Xyz: the whole text file as a data array
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, August 2018
%-------------------------------------------------------------------------

formatSpec = '%s%[^\n\r]';
fileID = fopen(Filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '',  'ReturnOnError', false);
dataArray{1} = strtrim(dataArray{1});
fclose(fileID);
Xyz = dataArray{:, 1};

Nelec = regexp((Xyz{1}),'\s','split'); Nelec = str2double(Nelec{1});
if (size(Xyz,1)-1)~=Nelec
    error('Number of electrodes in header does not match number of subsequent rows')
end
name = cell(Nelec,1);
x = nan(Nelec,1);
y = nan(Nelec,1);
z = nan(Nelec,1);
Count = 0;
for l = 2:size(Xyz,1)

    Parsed = regexp(Xyz(l,:),'\s','split');
    CurRow = Parsed{1}(~cellfun(@isempty,Parsed{1}));
    if numel(CurRow)~=4
        error('Expected x, y, z & electrode name but format does not match')
    else
        Count = Count+1;
        x(Count) = str2double(CurRow{1});
        y(Count) = str2double(CurRow{2});
        z(Count) = str2double(CurRow{3});
        name{Count} = CurRow{4};
    end
end

end
