function dual_merge_mrk( MrkFile1, MrkFile2 )
% dual_merge_mrk: merge markers in 2 .mrk files
%
% dual_merge_mrk( MrkFile1, MrkFile2 )
%
%  Inputs
% --------
% MrkFile1: first .mrk file 
% MrkFile2: second .mrk file 
%
%  Outputs
% ---------
% Both .mrk files with markers of the two files
% Original files will be backed up with suffix "_BKP"
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, October 2018
%-------------------------------------------------------------------------

% backup original files in case:
copyfile(MrkFile1,spm_file(MrkFile1,'suffix','_BKP'));
copyfile(MrkFile2,spm_file(MrkFile2,'suffix','_BKP'));

try
    [tS1, tE1, tL1] = read_mrk_file(MrkFile1);
    [tS2, tE2, tL2] = read_mrk_file(MrkFile2);
catch ME
    warning([ME.message,' Trying fallback marker reading function instead.'])
    [tS1, tE1, tL1] = read_mrk_Cartool(MrkFile1);
    [tS2, tE2, tL2] = read_mrk_Cartool(MrkFile2);
end

% if there are markers with identical labels in the two files, we need to
% be able to differentiate between the two (by adding a suffix e.g.):
CommonMarkers = intersect(unique(tL1),unique(tL2));
if ~isempty(CommonMarkers)
    for l = 1:length(CommonMarkers)
        tL1 = regexprep(tL1,CommonMarkers{l},['"',regexprep(CommonMarkers{l},'"',''),'_',spm_file(MrkFile1,'basename'),'"']);
        tL2 = regexprep(tL2,CommonMarkers{l},['"',regexprep(CommonMarkers{l},'"',''),'_',spm_file(MrkFile2,'basename'),'"']);
    end
end

% add markers where files were merged
FinalMarkerStart = [tS1;tS2];
FinalMarkerEnd = [tE1;tE2];
FinalMarkerLabel = [tL1;tL2];
[~,IdxSort] = sort(FinalMarkerStart);

write_mrk_file_Cartool(MrkFile1,FinalMarkerStart(IdxSort),FinalMarkerEnd(IdxSort),FinalMarkerLabel(IdxSort));
write_mrk_file_Cartool(MrkFile2,FinalMarkerStart(IdxSort),FinalMarkerEnd(IdxSort),FinalMarkerLabel(IdxSort));

end

