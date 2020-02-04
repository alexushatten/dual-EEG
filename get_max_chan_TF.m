function [ChanIdx,TFidx,MaxVal,ChanLabel,ChanCoor] = get_max_chan_TF(EEG,ChannelsInfo,Bipolar)
% GET_MAX_CHAN_TF: get channel and time frame with maximal value in EEG
% array, together with maximum value, and optionally, output channel label
% based on .els, .xyz or .spi file.
%
% [ChanIdx, TFidx, MaxVal] = get_max_chan_TF(EEG)
% [ChanIdx, TFidx, MaxVal, ChanLabel, ChanCoor] = ...
%                               get_max_chan_TF(EEG, ChannelsInfo)
% [ChanIdx, TFidx, MaxVal, ChanLabel, ChanCoor] = ...
%                               get_max_chan_TF(EEG, ChannelsInfo, Bipolar)
%
%  Inputs
% --------
% EEG:                      [channels x time] EEG traces
% ChannelsInfo (optional):  char, path to file embedding channels
%                           information, such as .els, .xyz, .spi.
%                           Alternatively, can be a cell array of strings
%                           with channel labels, whose size is equal to
%                           size(EEG,1).
% Bipolar (optional):       logical, whether the input EEG is bipolar (and
%                           the corresponding ChannelsInfo are not), in
%                           which case the ChannelsInfo will be converted
%                           to bipolar and the coordinates will be averaged
%                           by pairs (default = false).
%
%  Outputs
% ---------
% ChanIdx:   Index of channel with maximal value
% TFidx:     Index of time frame with maximal value
% MaxVal:    Value of maximum
% ChanLabel: Label of channel with maximal value, if ChannelsInfo was provided
% ChanCoor:  Coordinates of channel with maximal value, if ChannelsInfo was
%            provided
%
% See also GET_CHAN_COORDINATES
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, October 2019
%-------------------------------------------------------------------------

[ChanIdx,TFidx] = find(EEG == max(EEG(:)));
MaxVal = EEG(ChanIdx,TFidx);

if nargin < 2 && nargout > 3
    error('Output ChanLabel & ChanCoor requires ChannelsInfo')
end
if nargin < 3
    Bipolar = false;
end

if nargin > 1
    if isa(ChannelsInfo,'cell')
        if numel(ChannelsInfo) == size(EEG,1)
            ChanLabel = ChannelsInfo(ChanIdx);
        else
            error('Number of channels and labels do not match')
        end
    elseif isa(ChannelsInfo,'char')
        if exist(ChannelsInfo,'file')~=2
            error('File "%s" could not be found',ChannelsInfo)
        else
            if ~Bipolar
                % Make sure the number of channels in coordinates file and
                % input EEG match:
                [~,ChanLabel] = get_chan_coordinates(ChannelsInfo,1:size(EEG,1));
                if numel(ChanLabel)~=size(EEG,1)
                    error('Number of channels in input EEG and %s do not match',ChannelsInfo)
                else
                    [ChanCoor,ChanLabel] = get_chan_coordinates(ChannelsInfo,ChanIdx);
                end
            else
                if ~strcmpi(spm_file(ChannelsInfo,'ext'),'els')
                    error('Bipolar montage needs .els input file')
                else
                    [x,y,z,~,~,~,FullName] = read_els_file(ChannelsInfo);
                    [Bipoles,~,labelsB,labelsB2] = bipolar_montage(FullName,1);
                    if numel(labelsB) == size(EEG,1)
                        XYZ_bip = nan(size(labelsB,2),3);
                        for c = 1:size(Bipoles,1)
                            XYZ_bip(c,:) = mean([[x(Bipoles(c)),y(Bipoles(c)),z(Bipoles(c))];[x(Bipoles(c)+1),y(Bipoles(c)+1),z(Bipoles(c)+1)]]);
                        end
                        ChanCoor = XYZ_bip(ChanIdx);
                        ChanLabel = labelsB2(ChanIdx,:);
                    else
                        error('Number of generated labels and input channels do not match')
                    end
                end
            end
        end
    else
        error('ChannelsInfo provided but format is not recognized!')
    end
end

end

