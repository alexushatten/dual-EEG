function bipolar_scalp_AnyWaveMTG( Filename, Banana )
% bipolar_scalp_AnyWaveMTG: change AnyWave montage (.mtg) for scalp EEG to
% bipolar (double or triple banana)
%
% bipolar_scalp_AnyWaveMTG( Filename )
%
%  Inputs
% --------
% Filename: path to original .mtg file (XML-like file generated by AnyWave)
%           (must contain all channels to include in montage!)
% Banana: string, either 'double' or 'triple'
%
%  Outputs
% ---------
% Rewritten .mtg AnyWave file, original is kept but with ".bkp" suffix (if
% it does not exist already, otherwise it will simply overwrite the
% original file)
%
%-------------------------------------------------------------------------
% NB: modalities will be sorted alphabetically in the new .mtg file!
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% http://meg.univ-amu.fr/wiki/AnyWave
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, February 2019
%-------------------------------------------------------------------------

Montage = parseAnyWaveMTG(Filename);
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

% load('E:\code\MATLAB\MATLAB\dualEEG\hdchs.mat')
% MtgTxt = read_mtg_file('E:\code\MATLAB\MATLAB\AnyWave_display\For_making_template_montage\DoubleBanana_T7T8_lowT_right_left_OK.mtg');
hdchs = {'1';'F8';'3';'4';'F2';'6';'7';'8';'9';'AF8';'11';'AF4';'13';'14';'FCz';'16';'17';'FP2';'19';'20';'Fz';'22';'23';'FC1';'25';'FPz';'27';'28';'F1';'30';'31';'32';'33';'AF3';'35';'F3';'FP1';'38';'39';'40';'41';'FC3';'43';'C1';'45';'AF7';'F7';'F5';'FC5';'50';'51';'52';'53';'54';'55';'56';'57';'58';'C3';'60';'61';'FT7';'63';'C5';'65';'CP3';'FT9';'T9';'T7';'70';'71';'72';'73';'74';'75';'CP5';'77';'78';'CP1';'80';'81';'82';'83';'TP7';'85';'P5';'P3';'P1';'89';'CPz';'91';'92';'93';'TP9';'95';'P7';'PO7';'98';'99';'100';'Pz';'102';'103';'104';'105';'P9';'107';'108';'PO3';'110';'111';'112';'113';'114';'115';'O1';'117';'118';'POz';'120';'121';'122';'123';'124';'125';'Oz';'127';'128';'129';'130';'131';'132';'133';'134';'135';'136';'137';'138';'139';'PO4';'141';'P2';'CP2';'144';'145';'146';'147';'148';'149';'O2';'151';'152';'P4';'154';'155';'156';'157';'158';'159';'160';'PO8';'P6';'163';'CP4';'165';'166';'167';'168';'P10';'P8';'171';'CP6';'173';'174';'175';'176';'177';'178';'TP8';'180';'181';'182';'C4';'184';'C2';'186';'187';'188';'189';'TP10';'191';'192';'193';'C6';'195';'196';'197';'198';'199';'200';'201';'T8';'203';'204';'205';'FC4';'FC2';'208';'209';'T10';'FT8';'212';'FC6';'214';'215';'216';'217';'218';'FT10';'220';'221';'F6';'223';'F4';'225';'F10';'227';'228';'229';'230';'231';'232';'233';'234';'235';'236';'237';'238';'239';'240';'241';'242';'243';'244';'245';'246';'247';'248';'249';'250';'251';'F9';'253';'254';'255';'256';'Cz'};
MtgDouble = {'FP2','F10';'F10','T10';'T10','TP10';'TP10','P10';'P10','O2';'Fp2','F8';'F8','T8';'T8','P8';'P8','O2';'Fp2','F4';'F4','C4';'C4','P4';'P4','O2';'Fp1','F3';'F3','C3';'C3','P3';'P3','O1';'Fp1','F7';'F7','T7';'T7','P7';'P7','O1';'FP1','F9';'F9','T9';'T9','TP9';'TP9','P9';'P9','O1';'Fz','Cz';'Cz','Pz'};
MtgTriple = {'FP2','F10';'F10','T10';'T10','TP10';'TP10','P10';'P10','O2';'Fp2','F8';'F8','T8';'T8','P8';'P8','O2';'Fp2','F4';'F4','C4';'C4','P4';'P4','O2';'Fp1','F3';'F3','C3';'C3','P3';'P3','O1';'Fp1','F7';'F7','T7';'T7','P7';'P7','O1';'FP1','F9';'F9','T9';'T9','TP9';'TP9','P9';'P9','O1';'Fz','Cz';'Cz','Pz';'F10','F8';'F8','F4';'F4','Fz';'Fz','F3';'F3','F7';'F7','F9';'T10','T8';'T8','C4';'C4','Cz';'Cz','C3';'C3','T7';'T7','T9';'TP10','P8';'P8','P4';'P4','Pz';'Pz','P3';'P3','P7';'P7','TP9'};

switch Banana
    case 'double'
        Mtg = MtgDouble;
    case 'triple'
        Mtg = MtgTriple;
end

ElecGrp1 = match_vectors(strtrim(Mtg(:,1)),strtrim(hdchs),0); % case insensitive because "FP" in .xyz vs. "Fp" in .mtg
ElecGrp2 = match_vectors(strtrim(Mtg(:,2)),strtrim(hdchs),0); % usage of strtrim because some electrodes names in .mtg file have spaces appended

ScalpEEGchans = ~cellfun(@isempty,regexp(ChanType,'^EEG'));
IdxScalpEEGchans = find(ScalpEEGchans);
CheckMissing = match_vectors(hdchs,ChanName(ScalpEEGchans),1);
if isa(CheckMissing,'cell')
    warning('Some channel names are missing in AnyWave montage:')
    hdchs(cellfun(@isempty,CheckMissing))
end
AW2car2l = match_vectors(ChanName(ScalpEEGchans),hdchs,1);
if isa(AW2car2l,'cell')
    error('Unknown channel names in AnyWave montage!')
end

% ChanToKeep = ~cellfun(@isempty,match_vectors(ChanName,unique(hdchs([ElecGrp1;ElecGrp2])),1));
IdxChans1 = match_vectors(hdchs(ElecGrp1),ChanName(IdxScalpEEGchans),0);

%% We have to reconstruct the whole montage from scratch instead of just excluding what we don't need
% because the montage might be redundant (the same channel appears several times)

Modalities = unique(ChanType);
Counter = 1;
NewChanName = {};
NewChanType = {};
NewChanRef = {};
NewChanColor = {};
if ~NoFiltersSet
    NewChanFiltLP = [];
    NewChanFiltHP = [];
end
for m = 1:length(Modalities)
    if strcmp(Modalities(m),'EEG')
        NewChanName(Counter:Counter+length(ChanName(IdxScalpEEGchans(IdxChans1)))-1) ...
            = ChanName(IdxScalpEEGchans(IdxChans1));
        NewChanType(Counter:Counter+length(ChanName(IdxScalpEEGchans(IdxChans1)))-1) ...
            = ChanType(IdxScalpEEGchans(IdxChans1));
        % Make bipolar montage for scalp:
        NewChanRef(Counter:Counter+length(ChanName(IdxScalpEEGchans(IdxChans1)))-1) ...
            = hdchs(ElecGrp2);
        NewChanColor(Counter:Counter+length(ChanName(IdxScalpEEGchans(IdxChans1)))-1) ...
            = ChanColor(IdxScalpEEGchans(IdxChans1));
        if ~NoFiltersSet
            NewChanFiltLP(Counter:Counter+length(ChanName(IdxScalpEEGchans(IdxChans1)))-1) ...
                = ChanFiltLP(IdxScalpEEGchans(IdxChans1));
            NewChanFiltHP(Counter:Counter+length(ChanName(IdxScalpEEGchans(IdxChans1)))-1) ...
                = ChanFiltHP(IdxScalpEEGchans(IdxChans1));
        end
    else
        NewChanName(Counter:Counter+sum(strcmp(ChanType,Modalities(m)))-1) ...
            = ChanName(strcmp(ChanType,Modalities(m)));
        NewChanType(Counter:Counter+sum(strcmp(ChanType,Modalities(m)))-1) ...
            = ChanType(strcmp(ChanType,Modalities(m)));
        NewChanRef(Counter:Counter+sum(strcmp(ChanType,Modalities(m)))-1) ...
            = ChanRef(strcmp(ChanType,Modalities(m)));
        NewChanColor(Counter:Counter+sum(strcmp(ChanType,Modalities(m)))-1) ...
            = ChanColor(strcmp(ChanType,Modalities(m)));
        if ~NoFiltersSet
            NewChanFiltLP(Counter:Counter+sum(strcmp(ChanType,Modalities(m)))-1) ...
                = ChanFiltLP(strcmp(ChanType,Modalities(m)));
            NewChanFiltHP(Counter:Counter+sum(strcmp(ChanType,Modalities(m)))-1) ...
                = ChanFiltHP(strcmp(ChanType,Modalities(m)));
        end
    end
    Counter = Counter + sum(strcmp(NewChanType,Modalities(m)));
end

Montage.ChanName = NewChanName';
Montage.ChanType = NewChanType';
Montage.ChanRef = NewChanRef';
Montage.ChanColor = NewChanColor';
if ~NoFiltersSet
    Montage.ChanFiltLP = NewChanFiltLP';
    Montage.ChanFiltHP = NewChanFiltHP';
end

% Backup original .mtg
if exist(spm_file(Filename,'suffix','.bkp'),'file')~=2
    copyfile(Filename,spm_file(Filename,'suffix','.bkp'));
end

% Write new file
writeAnyWaveMTG(Filename,Montage);

end
