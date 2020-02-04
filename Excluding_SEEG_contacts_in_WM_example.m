
% In Cartool, change orientation of T1 that is aligned to CT and SEEG electrodes and
% set origin to 0,0,0

% Perform segmentation of T1 in SPM (bias-corrected T1 can also be
% requested to check results manually and tissue segmentation results!)

% Select .els files that contains all SEEG contacts:
[x,y,z,name,ClusterName,Nclus,FullName,Els] = read_els_file(spm_select);
% VmT1 = spm_vol(spm_select); % select bias-corrected T1
Vc2T1 = spm_vol(spm_select); % select tissue class #2 (WM)
% DmT1 = spm_read_vols(VmT1); % read the data
Dc2T1 = spm_read_vols(Vc2T1); % read the data
ElC2 = nan(size(x));
for c = 1:length(x)
    ElC2(c) = Dc2T1(round(x(c))+1,round(y(c))+1,round(z(c))+1);
end
% ElT1 = nan(size(x));
% for c = 1:length(x)
%     ElT1(c) = DmT1(round(x(c))+1,round(y(c))+1,round(z(c))+1);
% end
% To be removed:
sum(ElC2>((256*0.99)-1)) % this can be adapted based on the results (one can try 0.75, 0.9, 0.95...)
figure; hist(ElC2,100)

% Will be excluded:
FullName(ElC2>((256*0.99)-1))'

% To be kept:
sum(~(ElC2>((256*0.99)-1)))

% Will be kept:
FullName(~(ElC2>((256*0.99)-1)))'

% However, we don't want to forget about pairs for which one contact is in
% GM and the other is in WM, so we will dilate this a little...
ToKeepTemp = (conv([~(ElC2>((256*0.99)-1))]'+0,[1 1 1]')>0);
ToKeep = ToKeepTemp(2:end-1);

ToRemove = ~ToKeep;
% 39 contacts to exclude for SPES (CCEP) !

[Bipoles,ElecMatch,labelsB,labelsB2] = bipolar_montage(FullName(ToKeep),1);



