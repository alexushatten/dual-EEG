function colors_AnyWave_mrk( Filename, ToInclude, ToIgnore )
% colors_AnyWaveMTG_mrk: change default color of markers in AnyWave
%
% colors_AnyWaveMTG_mrk( Filename, ToInclude, ToIgnore )
%
%  Inputs
% --------
% Filename:  path to AnyWave .mrk file
% ToInclude: markers for which color should be added based on their unique type,
%            regex allowed, e.g. starting with "manual_" and ending with
%            "pos" would be '^manual_.*pos$', default = '.*' (all)
% ToIgnore:  markers that should be ignored for making the color lookup table
%
%  Outputs
% ---------
% Rewritten .mrk AnyWave file, original is kept but with ".bkp" suffix (if
% it does not exist already, otherwise it will simply overwrite the
% original file)
%
%-------------------------------------------------------------------------
% http://meg.univ-amu.fr/wiki/AnyWave
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, February 2019
%-------------------------------------------------------------------------

if nargin < 2
    ToInclude = '.*';
end
if nargin < 3
    ToIgnore = '';
end

[AWcolorCode, AWcolorName, AWlightcolorNice] = get_color_names_and_codes;

% Exclude black because it is default color:
AWlightcolorNice = AWlightcolorNice(cellfun(@isempty,regexp(AWlightcolorNice,'black')));

AWcolorCodeNice = AWcolorCode(~cellfun(@isempty,match_vectors(AWcolorName,AWlightcolorNice,0)));

[Label,Code,Onset,Duration,Comment,Color] = read_mrk_AW(Filename);

if sum(cellfun(@isempty,Color))<length(Label)
    OverwritePreExistingColors = questdlg('There are already colors set in this .mrk file, are you sure you want to overwrite?','Warning: overwrite pre-existing colors?','Yes','No','No');
    if strcmp(OverwritePreExistingColors,'No')
        warning('Leaving the existing .mrk file untouched and exiting...')
        return;
    end
end

% Make groups based on ToInclude:
InterestingLabels = Label(~cellfun(@isempty,regexp(Label,ToInclude)));

% Filter out labels to ignore:
if ~isempty(ToIgnore)
    InterestingLabels = InterestingLabels(cellfun(@isempty,regexp(InterestingLabels,ToIgnore)));
end

% Lookup table:
LUT = [unique(InterestingLabels),AWcolorCodeNice(1:length(unique(InterestingLabels)))];

for l = 1:length(Label)
    LUTidx = match_vectors(Label(l),LUT(:,1),1);
    if ~isa(LUTidx,'cell')
        Color(l) = LUT(LUTidx,2);
    end
end

if exist(spm_file(Filename,'suffix','.bkp'),'file')~=2
    copyfile(Filename,spm_file(Filename,'suffix','.bkp'));
end

write_mrk_file_AnyWave(Filename,Label,Onset,Duration,Comment,Color,Code)

end

function [AWcolorCode, AWcolorName, AWcolorNice] = get_color_names_and_codes

% The following set of colors should be different from dark colors as
% defined in colors_groups_SEEG_AnyWaveMTG.m (because these are best for
% channels, but for markers light colors are likely more appropriate)
%
% Nevertheless, we should avoid magenta, fuchsia, and everything close to
% pink because this is the default value for markers with duration > 0...
AWcolorNice = {'aqua'
    'darksalmon'
    'aquamarine'
    'khaki'
    'chartreuse'
    'deepskyblue'
    'springgreen'
    'gold'
    'lightgreen'
    'darkseagreen'
    'darkturquoise'
    'lime'
    'mediumaquamarine'
    'moccasin'
    'palegreen'
    'navajowhite'
    'mediumspringgreen'
    'thistle'
    'tomato'
    'paleturquoise'
    'skyblue'
    'greenyellow'
    'cyan'
    'peru'
    'powderblue'
    'orange'
    'lightcoral'
    'lawngreen'
    'bisque'
};

AWcolorCode = {'#f0f8ff'
'#faebd7'
'#00ffff'
'#7fffd4'
'#f0ffff'
'#f5f5dc'
'#ffe4c4'
'#000000'
'#ffebcd'
'#0000ff'
'#8a2be2'
'#a52a2a'
'#deb887'
'#5f9ea0'
'#7fff00'
'#d2691e'
'#ff7f50'
'#6495ed'
'#fff8dc'
'#dc143c'
'#00ffff'
'#00008b'
'#008b8b'
'#b8860b'
'#a9a9a9'
'#006400'
'#a9a9a9'
'#bdb76b'
'#8b008b'
'#556b2f'
'#ff8c00'
'#9932cc'
'#8b0000'
'#e9967a'
'#8fbc8f'
'#483d8b'
'#2f4f4f'
'#2f4f4f'
'#00ced1'
'#9400d3'
'#ff1493'
'#00bfff'
'#696969'
'#696969'
'#1e90ff'
'#b22222'
'#fffaf0'
'#228b22'
'#ff00ff'
'#dcdcdc'
'#f8f8ff'
'#ffd700'
'#daa520'
'#808080'
'#008000'
'#adff2f'
'#808080'
'#f0fff0'
'#ff69b4'
'#cd5c5c'
'#4b0082'
'#fffff0'
'#f0e68c'
'#e6e6fa'
'#fff0f5'
'#7cfc00'
'#fffacd'
'#add8e6'
'#f08080'
'#e0ffff'
'#fafad2'
'#d3d3d3'
'#90ee90'
'#d3d3d3'
'#ffb6c1'
'#ffa07a'
'#20b2aa'
'#87cefa'
'#778899'
'#778899'
'#b0c4de'
'#ffffe0'
'#00ff00'
'#32cd32'
'#faf0e6'
'#ff00ff'
'#800000'
'#66cdaa'
'#0000cd'
'#ba55d3'
'#9370db'
'#3cb371'
'#7b68ee'
'#00fa9a'
'#48d1cc'
'#c71585'
'#191970'
'#f5fffa'
'#ffe4e1'
'#ffe4b5'
'#ffdead'
'#000080'
'#fdf5e6'
'#808000'
'#6b8e23'
'#ffa500'
'#ff4500'
'#da70d6'
'#eee8aa'
'#98fb98'
'#afeeee'
'#db7093'
'#ffefd5'
'#ffdab9'
'#cd853f'
'#ffc0cb'
'#dda0dd'
'#b0e0e6'
'#800080'
'#ff0000'
'#bc8f8f'
'#4169e1'
'#8b4513'
'#fa8072'
'#f4a460'
'#2e8b57'
'#fff5ee'
'#a0522d'
'#c0c0c0'
'#87ceeb'
'#6a5acd'
'#708090'
'#708090'
'#fffafa'
'#00ff7f'
'#4682b4'
'#d2b48c'
'#008080'
'#d8bfd8'
'#ff6347'
'#000000'
'#40e0d0'
'#ee82ee'
'#f5deb3'
'#ffffff'
'#f5f5f5'
'#ffff00'
'#9acd32'};

AWcolorName = {'aliceblue';'antiquewhite';'aqua';'aquamarine';'azure';'beige';'bisque';'black';'blanchedalmond';'blue';'blueviolet';'brown';'burlywood';'cadetblue';'chartreuse';'chocolate';'coral';'cornflowerblue';'cornsilk';'crimson';'cyan';'darkblue';'darkcyan';'darkgoldenrod';'darkgray';'darkgreen';'darkgrey';'darkkhaki';'darkmagenta';'darkolivegreen';'darkorange';'darkorchid';'darkred';'darksalmon';'darkseagreen';'darkslateblue';'darkslategray';'darkslategrey';'darkturquoise';'darkviolet';'deeppink';'deepskyblue';'dimgray';'dimgrey';'dodgerblue';'firebrick';'floralwhite';'forestgreen';'fuchsia';'gainsboro';'ghostwhite';'gold';'goldenrod';'gray';'green';'greenyellow';'grey';'honeydew';'hotpink';'indianred';'indigo';'ivory';'khaki';'lavender';'lavenderblush';'lawngreen';'lemonchiffon';'lightblue';'lightcoral';'lightcyan';'lightgoldenrodyellow';'lightgray';'lightgreen';'lightgrey';'lightpink';'lightsalmon';'lightseagreen';'lightskyblue';'lightslategray';'lightslategrey';'lightsteelblue';'lightyellow';'lime';'limegreen';'linen';'magenta';'maroon';'mediumaquamarine';'mediumblue';'mediumorchid';'mediumpurple';'mediumseagreen';'mediumslateblue';'mediumspringgreen';'mediumturquoise';'mediumvioletred';'midnightblue';'mintcream';'mistyrose';'moccasin';'navajowhite';'navy';'oldlace';'olive';'olivedrab';'orange';'orangered';'orchid';'palegoldenrod';'palegreen';'paleturquoise';'palevioletred';'papayawhip';'peachpuff';'peru';'pink';'plum';'powderblue';'purple';'red';'rosybrown';'royalblue';'saddlebrown';'salmon';'sandybrown';'seagreen';'seashell';'sienna';'silver';'skyblue';'slateblue';'slategray';'slategrey';'snow';'springgreen';'steelblue';'tan';'teal';'thistle';'tomato';'transparent';'turquoise';'violet';'wheat';'white';'whitesmoke';'yellow';'yellowgreen'};

end
