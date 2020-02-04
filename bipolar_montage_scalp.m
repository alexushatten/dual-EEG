function [bip_labels,bip_labels_m,ElecGrp1,ElecGrp2,bip_EEG] = bipolar_montage_scalp(Banana,EEG)
% Make bipolar montage from electrode names and return EEG trace with
% bipolar montage if entered as input (for SEEG, based on channel labels)
% 
% Usage:
%-------------------------------------------------------------------------
% [bip_labels,bip_labels_m,ElecGrp1,ElecGrp2] = bipolar_montage_scalp(Banana)
%
% [bip_labels,bip_labels_m,ElecGrp1,ElecGrp2,bip_EEG] = bipolar_montage_scalp(Banana,EEG)
% 
% Inputs:
%-------------------------------------------------------------------------
% Banana: character, 'double' or 'triple' based on the montage desired
%
% EEG: (optional) [channel x time] EEG traces to recalculate based on bipolar
% montage to be created
%
% Outputs:
%-------------------------------------------------------------------------
% bip_labels: [n x 2] cell array with labels of channel pairs
%
% bip_labels_m: [n x 1] cell array with pairs of channels separated by "-"
%
% ElecGrp1: left hand side of electrodes included in montage (to be substracted from)
%
% ElecGrp2: right hand side of electrodes included in montage (substracted to ElecGrp1)
% 
% bip_EEG: [channel pairs x time] EEG traces recalculated based on bipolar
% montage
%-------------------------------------------------------------------------
% NB: for high density (257 channels) EEG only!
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, April 2019
%-------------------------------------------------------------------------

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

bip_labels = [hdchs(ElecGrp1),hdchs(ElecGrp2)];

clear bip_labels_m
for c = 1:size(bip_labels,1)
    bip_labels_m{c} = [bip_labels{c,1},'-',bip_labels{c,2}]; %#ok<AGROW>
end

if nargin>1
    bip_EEG = nan(length(ElecGrp1),size(EEG,2));
    for c = 1:size(bip_EEG,1)
        bip_EEG(c,:) = EEG(ElecGrp1(c),:)-EEG(ElecGrp2(c),:);
    end
end

end
