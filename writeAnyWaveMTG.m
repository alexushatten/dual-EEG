function writeAnyWaveMTG( Filename, Montage)
% writeAnyWaveMTG: writes an AnyWave montage (.mtg) file
%
% Struct = parseAnyWaveMTG( Filename )
%
%  Inputs
% --------
% Filename   : filename to write
% Montage    : structure with fields:
%   ChanName   : [n x 1] cell of strings with channel names
%   ChanType   : [n x 1] cell of strings with channel types
%   ChanRef    : [n x 1] cell of strings with reference for each channel
%   ChanColor  : [n x 1] cell of strings with channel colors
%   ChanFiltLP (optional) : [n x 1] array with low-pass filter for each channel
%   ChanFiltHP (optional) : [n x 1] array with high-pass filter for each channel
%
%  Output
% --------
% AnyWave's XML-like montage (.mtg) file based on input parameters
%
%-------------------------------------------------------------------------
% http://meg.univ-amu.fr/wiki/AnyWave
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, February 2019
%-------------------------------------------------------------------------

ChanName = Montage.ChanName;
ChanType = Montage.ChanType;
ChanRef = Montage.ChanRef;
ChanColor = Montage.ChanColor;
if isfield(Montage,'ChanFiltLP') && isfield(Montage,'ChanFiltHP')
    ChanFiltLP = Montage.ChanFiltLP;
    ChanFiltHP = Montage.ChanFiltHP;
    NoFiltersSet = false;
else
    NoFiltersSet = true;
end

Head{1} = '<!DOCTYPE AnyWaveMontage>';
Head{2} = '<Montage>';
Head = Head';

% ====== TEMPLATE ======
% '<Channel name="CHANNEL1">'
% '<type>EEG</type>'
% '<reference>AVG</reference>'
% '<color>black</color>'
% '<filters lowPass="70" highPass="1"/>'
% '</Channel>'
% ======================

Counter = 1;
for c = 1:length(ChanName)
    Body{Counter} = ['   <Channel name="',ChanName{c},'">']; %#ok<*AGROW,*NASGU>
    Counter = Counter + 1;
    Body{Counter} = ['      <type>',ChanType{c},'</type>'];
    Counter = Counter + 1;
    Body{Counter} = ['      <reference>',ChanRef{c},'</reference>'];
    Counter = Counter + 1;
    Body{Counter} = ['      <color>',ChanColor{c},'</color>'];
    if ~NoFiltersSet
        Counter = Counter + 1;
        Body{Counter} = ['      <filters highPass="',num2str(ChanFiltHP(c)),'" lowPass="',num2str(ChanFiltLP(c)),'"/>'];
    end
    Counter = Counter + 1;
    Body{Counter} = '   </Channel>';
    Counter = Counter + 1;
end
Body = Body';

Tail = {'</Montage>'};
Contents = [Head;Body;Tail];

fileID = fopen(Filename,'w');
for l = 1:size(Contents,1)
    fprintf(fileID,'%s\n',Contents{l});
end
fclose(fileID);

end

