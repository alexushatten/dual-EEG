
AllPathsSEF = {'paths_to_all_SEF_files_for_which_ADES_has_to_be_fixed'};

AllPathsMRK = spm_file(AllPathsSEF,'ext','sef.mrk');

load('E:\code\MATLAB\MY_CODE\matlab\dualEEG\hdchs.mat')

for f = 1:2:(length(AllPathsSEF)-1)
    fprintf('Doing file %d/%d...\n',f,length(AllPathsSEF))
    
    load(spm_file(AllPathsSEF{f},'ext','mat'));
    load(spm_file(AllPathsSEF{f+1},'ext','mat'));
    try
        [hd_T1_aligned_ALL,hd_T2_aligned_ALL,hd_L_aligned_ALL] = read_mrk_Cartool(AllPathsMRK{f});
    catch
        [hd_T1_aligned_ALL,hd_T2_aligned_ALL,hd_L_aligned_ALL] = read_mrk_Cartool(spm_file(AllPathsMRK{f},'suffix','_BKP_before_marking'));
    end
    try
        [ic_T1_aligned_ALL,ic_T2_aligned_ALL,ic_L_aligned_ALL] = read_mrk_Cartool(AllPathsMRK{f+1});
    catch
        [ic_T1_aligned_ALL,ic_T2_aligned_ALL,ic_L_aligned_ALL] = read_mrk_Cartool(spm_file(AllPathsMRK{f+1},'suffix','_BKP_before_marking'));
    end
    mat2ades([ic_EEGialigned;hd_EEGaligned],...
        spm_file(regexprep(AllPathsSEF{f},'hdEEG','dual'),'ext',''),...
        hd_Hdr.samplingfreq,[ic_Hdr.channelnames';hdchs],...
        [cellstr(repmat('SEEG',length(ic_Hdr.channelnames),1));...
        cellstr(repmat('EEG',length(hdchs),1))]);
    
    Ms = [ic_T1_aligned_ALL;hd_T1_aligned_ALL];
    Me = [ic_T2_aligned_ALL;hd_T2_aligned_ALL];
    ML = [ic_L_aligned_ALL;hd_L_aligned_ALL];
    [~,SortIdx] = sort(Ms);
    write_mrk_file_AnyWave([spm_file(regexprep(AllPathsSEF{f},'hdEEG','dual'),'ext',''),'.ades.mrk'],...
        ML(SortIdx),Ms(SortIdx)/hd_Hdr.samplingfreq,(Me(SortIdx)-Ms(SortIdx))/hd_Hdr.samplingfreq);
    
end

