function comment_AnyWave_mrk( Filename, MrkLabels, ChanLabels )
% colors_AnyWaveMTG_mrk: add comment to markers in AnyWave based on
% channel labels
%
% colors_AnyWaveMTG_mrk( Filename, ToInclude, ToIgnore )
%
%  Inputs
% --------
% Filename:               path to AnyWave .mrk file
% MrkLabels:              cell array of strings with marker types for which
%                         to add comments, e.g. {'HAG12_pos', 'HPG23_neg'}
% ChanLabels [optional] : cell array of strings with all channels
%                         implicated for each marker types, separated with
%                         commas when there are more than one, e.g.:
%                         {'HAG1-HAG2', 'HPG1-HPG2,HPG2-HPG3'}
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
% Renaud Marquis @ FBMlab, March 2019
%-------------------------------------------------------------------------

[Label,Code,Onset,Duration,Comment,Color] = read_mrk_AW(Filename);
ChanLabels = ChanLabels(:);

for m = 1:length(Label)
    TempMatch = match_vectors(Label(m),MrkLabels,1);
    if ~isa(TempMatch,'cell')
        Comment{m} = [ChanLabels{TempMatch},','];
    end
end

if exist(spm_file(Filename,'suffix','.bkp'),'file')~=2
    copyfile(Filename,spm_file(Filename,'suffix','.bkp'));
end

write_mrk_file_AnyWave(Filename,Label,Onset,Duration,Comment,Color,Code)

end

% % Alternative: parse .ades file and/or .mtg file if it exists and match
% MrkLabels to it, but this is too complex because even if we feed in the
% relevant markers only (and not all of them), not all channels have an
% associated marker, and some channels have arbitrary numbers (not only
% SEEG with electrodes going up to >9 but also scalp EEG, ECG, MRK+/-
% channels...), so it would be very complex to make it work and not
% necessarily very practical, so likely better to simply feed in the
% ChanLabels as well, because the hard work is to add the Comment
% ("Target") to each marker, which will still be done using the code above.

% First check if there exist a .mtg file, because in such case the Targets
% should be defined based on that!

% [RootPath,ADESfile] = fileparts(Filename);
% 
% if ~(strcmp(regexprep(Filename,'.mrk',''),fullfile(RootPath,ADESfile)) || strcmp(Filename(1:end-4),fullfile(RootPath,ADESfile)))
%     warning('File name / extension of metadata format might be incorrect...')
% end
% 
% if nargin < 3
%     [Channels, Modality, Nsamples, FreqSamp] = read_ades_AW(fullfile(RootPath,ADESfile));
% end
% 
% % ===== Parser of labels based on regular expressions =====
% % special characters:
% for m = 1:length(MrkLabels)
%     ParsingStep1 = regexp(MrkLabels{m},'[^a-zA-Z0-9]','split');
%     Numbers = {};
%     for p = 1:length(ParsingStep1)
%         TempNumbers = regexp(ParsingStep1{p},'[0-9]','match');
%         if ~isempty(TempNumbers)
%             Numbers(end+1:end+length(TempNumbers)) = TempNumbers; 
%         end
%     end
%     Parts = regexprep(ParsingStep1,'[0-9]','');
%     Parts(end+1:end+length(Numbers)) = Numbers; 
%     MetaParts{m} = Parts; %#ok<AGROW,NASGU>
% end
% 
% for m = 1:length(MrkLabels)
%     
% end

% Exclude black because it is default color:
