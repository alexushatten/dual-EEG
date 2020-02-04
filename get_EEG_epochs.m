function Epochs = get_EEG_epochs(EEG,TimeFrames,Pre,Post)
% get_EEG_epochs: extract epochs from EEG traces
%
% get_EEG_epochs( EEG, TimeFrames )
%
%  Inputs
% --------
% EEG : [channels x time] EEG traces
% TimeFrames: [n x 1] time frames of markers
% Pre: [1 x 1] number of time frames to keep before each marker
% Post: [1 x 1] number of time frames to keep after each marker
%
%  Outputs
% ---------
% Epochs: [channels x [Pre+Post+1] x epochs] epochs from EEG traces
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, September 2019
%-------------------------------------------------------------------------

% Discard epochs too close to BOF or EOF
if any((TimeFrames-Pre)<1)
    warning('Some epochs were too close to beginning of EEG are were discarded!')
end
TimeFrames((TimeFrames-Pre)<1)=[];
if any((TimeFrames+Post)>size(EEG,2))
    warning('Some epochs were too close to end of EEG are were discarded!')
end
TimeFrames((TimeFrames+Post)>size(EEG,2))=[];

Epochs = nan(size(EEG,1),Pre+Post+1,length(TimeFrames));
for ep = 1:length(TimeFrames)
    Epochs(:,:,ep) = EEG(:,(TimeFrames(ep)-Pre):(TimeFrames(ep)+Post));
end

end

