function Epochs = markers2epochs( MarkerStart, MarkerEnd, EpochDuration, Method )
% markers2epochs: create subscript indices vector of epochs based on
% markers start and end vectors given an epoch duration
%
% Epochs = markers2epochs( MarkerStart, MarkerEnd, EpochDuration, Method )
%
%  Inputs
% --------
% MarkerStart:      [n x 1], subscript indices of markers's start time
% MarkerEnd:        [n x 1], subscript indices of markers's start end
% EpochDuration:    [1 x n], subscript indices of requested epochs around
%                   markers, e.g. -9:1:10
% Method:           string, what to do when for a given marker, MarkerStart
%                   and MarkerEnd are different (i.e. duration > 0):
%                   - consider marker start and end and add EpochDuration
%                   around it ('augment')
%                   - take the center of the pre-existing epoch (i.e. the timing
%                   between start and end) and add EpochDuration around it ('center')
%                   - consider marker start and ignore end ('start')
%                   - consider marker end and ignore start ('end')
%
%  Outputs
% ---------
% Epochs:           subscript indices of epochs around markers with duration
%                   length(EpochDuration), except for markers close to first
%                   sample and if Method = 'augment'
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, November 2018
%-------------------------------------------------------------------------

if strcmpi(Method,'center')
    Marker = round((MarkerStart+MarkerEnd)/2);
    Epochs = unique(vect(bsxfun(@plus,Marker,EpochDuration)));
elseif strcmpi(Method,'start')
    Epochs = unique(vect(bsxfun(@plus,MarkerStart,EpochDuration)));
elseif strcmpi(Method,'end')
    Epochs = unique(vect(bsxfun(@plus,MarkerEnd,EpochDuration)));
elseif strcmpi(Method,'augment')
    Epochs1 = sort(vect(bsxfun(@plus,MarkerStart,EpochDuration)));
    Epochs2 = sort(vect(bsxfun(@plus,MarkerEnd,EpochDuration)));
    Epochs = unique([Epochs1;Epochs2]);
end

% get rid of timings below 1:
Epochs(Epochs<1)=[];

end

