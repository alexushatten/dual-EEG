function [Coordinates,Label] = get_chan_coordinates(ChannelsInfo,ChannelsIndices,Avg)
% GET_CHAN_COORDINATES: get channel (sensor or solution point) coordinates
% and label based on electrode setup (.els), electrode coordinates (.xyz)
% or solution points (.spi) files
%
% [Coordinates, Label] = get_chan_coordinates( ChannelsInfo, ChannelsIndices,Avg )
%
%  Inputs
% --------
% ChannelsInfo: char, path to .els, .xyz or .spi file from Cartool
% ChannelsIndices: integer, index or vector of indices of channels
% Avg (optional): whether to average coordinates when length of
%                 ChannelsIndices > 1. Useful for bipolar montage of SEEG data
%                 (default = false).
%
%  Outputs
% ---------
% Coordinates: [1 x n] coordinates of channel, average coordinates if multiple
%              indices were requested and if Avg = true, otherwise [c x n]
%              where c is the number of channel indices and n is the number
%              of dimensions of coordinates in space (usually [c x 3])
% Label:       cell array of characters, label of channel index, or
%              concatenated labels of channels if multiple channel indices
%              and if Avg = true, using "-" used as separator, like e.g. :
%                   {'channel1-channel2'}
%              If Avg = false, {c x 1} cell array of chars.
%
% See also GET_MAX_CHAN_TF
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, October 2019
%-------------------------------------------------------------------------

if nargin < 3
    Avg = false;
end

Ext = spm_file(ChannelsInfo,'ext');
switch lower(Ext)
    case 'els'
        [x,y,z,~,~,~,FullName] = read_els_file(ChannelsInfo);
        XYZ = [x,y,z];
        Names = FullName;
    case 'xyz'
        [x,y,z,name] = read_xyz_file(ChannelsInfo);
        XYZ = [x,y,z];
        Names = name;
    case 'spi'
        SPs = read_spi_Cartool(ChannelsInfo);
        XYZ = [SPs.x,SPs.y,SPs.z];
        Names = SPs.labels;
end

if length(ChannelsIndices)>1 && Avg
    Coordinates = mean(XYZ(ChannelsIndices,:));
    Label = [];
    Count = 1;
    for c = 1:length(ChannelsIndices)
        L = length(Names(ChannelsIndices(c)));
        Label(Count:Count+L) = Names(ChannelsIndices(c));
        Label(end:end+1) = '-';
        Count = Count+L+1;
    end
    Label = cellstr(Label);
else
    Coordinates = XYZ(ChannelsIndices,:);
    Label = Names(ChannelsIndices);
end

end

