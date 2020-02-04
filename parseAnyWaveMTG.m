function Montage = parseAnyWaveMTG( Filename )
% parseAnyWaveMTG: parse an AnyWave montage (.mtg) file
%
% Montage = parseAnyWaveMTG( Filename )
%
%  Inputs
% --------
% Filename: path to .mtg file (XML-like file generated by AnyWave)
%
%  Outputs
% ---------
% Montage    : structure with fields:
%   ChanName   : [n x 1] cell of strings with channel names
%   ChanType   : [n x 1] cell of strings with channel types
%   ChanRef    : [n x 1] cell of strings with reference for each channel
%   ChanColor  : [n x 1] cell of strings with channel colors
%   ChanFiltLP (optional) : [n x 1] array with low-pass filter for each channel
%   ChanFiltHP (optional) : [n x 1] array with high-pass filter for each channel
%
%-------------------------------------------------------------------------
% http://meg.univ-amu.fr/wiki/AnyWave
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, February 2019
%-------------------------------------------------------------------------

Mtg = read_vrb_file(Filename);

%% Channel name
ChanNameString = '<Channel name=';
ChanNamePos = find(~cellfun(@isempty,regexp(Mtg,ChanNameString)));
Start = length(ChanNameString)+2;
End = cell2mat(regexp(Mtg(ChanNamePos),'">'))-1;
for c = 1:length(ChanNamePos)
    Temp = Mtg{ChanNamePos(c)};
    ChanName{c} = Temp(Start:End(c)); %#ok<*AGROW>
end
ChanName = ChanName';

%% Channel type
ChanTypeString = '<type>';
ChanTypePos = find(~cellfun(@isempty,regexp(Mtg,ChanTypeString)));
Start = length(ChanTypeString)+1;
End = cell2mat(regexp(Mtg(ChanTypePos),'</'))-1;
for c = 1:length(ChanTypePos)
    Temp = Mtg{ChanTypePos(c)};
    ChanType{c} = Temp(Start:End(c));
end
ChanType = ChanType';

%% Channel reference
ChanRefString = '<reference>';
ChanRefPos = find(~cellfun(@isempty,regexp(Mtg,ChanRefString)));
Start = length(ChanRefString)+1;
End = cell2mat(regexp(Mtg(ChanRefPos),'</'))-1;
for c = 1:length(ChanRefPos)
    Temp = Mtg{ChanRefPos(c)};
    ChanRef{c} = Temp(Start:End(c));
end
ChanRef = ChanRef';

%% Channel color
ChanColorString = '<color>';
ChanColorPos = find(~cellfun(@isempty,regexp(Mtg,ChanColorString)));
Start = length(ChanColorString)+1;
End = cell2mat(regexp(Mtg(ChanColorPos),'</'))-1;
for c = 1:length(ChanColorPos)
    Temp = Mtg{ChanColorPos(c)};
    ChanColor{c} = Temp(Start:End(c));
end
ChanColor = ChanColor';

%% Channel filters
ChanFiltString = '<filters lowPass="';
ChanFiltPos = find(~cellfun(@isempty,regexp(Mtg,ChanFiltString)));
ChanFiltParts = regexp(Mtg(ChanFiltPos),'"','split');
if ~isempty(ChanFiltParts)
    for c = 1:length(ChanFiltPos)
        Temp = ChanFiltParts{c};
        ChanFiltLP(c) = str2double(Temp(2));
        ChanFiltHP(c) = str2double(Temp(4));
    end
    ChanFiltLP = ChanFiltLP';
    ChanFiltHP = ChanFiltHP';
end

Montage.ChanName = ChanName;
Montage.ChanType = ChanType;
Montage.ChanRef = ChanRef;
Montage.ChanColor = ChanColor;
if ~isempty(ChanFiltParts)
    Montage.ChanFiltLP = ChanFiltLP;
    Montage.ChanFiltHP = ChanFiltHP;
end

end
