function write_mrk_file_AnyWave(Filename,Label,Onset,Duration,Comment,Color,Code)
% write_mrk_AW: writes AnyWave .mrk file
%
% write_mrk_file_AnyWave(Filename,Label,Onset)
% write_mrk_file_AnyWave(Filename,Label,Onset,Duration)
% write_mrk_file_AnyWave(Filename,Label,Onset,Duration,Comment)
% write_mrk_file_AnyWave(Filename,Label,Onset,Duration,Comment,Color)
% write_mrk_file_AnyWave(Filename,Label,Onset,Duration,Comment,Color,Code)
%
%  Inputs
% --------
% Filename:             path to AnyWave .mrk file, do not forget that
%                       filename/extension should end with .ades.mrk !
% Label:                cell array of strings with labels
% Onset:                numeric array with onsets in seconds
% Duration [optional]:  numeric array with durations in seconds, default = 0
% Comment [optional]:   cell of strings with comments, default = ''
% Color [optional]:     cell of strings with color code, default = '', see
%                       colors_AnyWave_mrk.m internal function for more details
% Code [optional]:      array with numeric code, default = -1
%
%  Outputs
% ---------
% AnyWave .mrk file @ path Filename
%
% http://meg.univ-amu.fr/wiki/AnyWave:ADES#The_marker_file_.28.mrk.29
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, October 2018
%                          March 2019: rewritten with more options
%-------------------------------------------------------------------------

if nargin < 4
    Duration = zeros(length(Label),1);
%     FlagNoDuration = true;
end
if nargin < 5
    Comment = cell(length(Label),1);
%     FlagNoComment = true;
end   
if nargin < 6
    Color = cell(length(Label),1);
%     FlagNoColor = true;
end
if nargin < 7
    Code = -ones(length(Label),1);
%     FlagNoCode = true;
end

fileID = fopen(Filename,'w');
fprintf(fileID,'%3s\r\n','// AnyWave Marker File');
for n = 1:size(Label,1)
    formatSpec = '%s\t%.0f\t%.5f\t%.5f\t%s\t%s\r\n';
    fprintf(fileID,formatSpec,Label{n},Code(n),Onset(n),Duration(n),Comment{n},Color{n});
end
fclose(fileID);

end
