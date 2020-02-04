
[~,~,mgrid_fp] = my_recursive_listfiles('E:\FS_subjects_DONE','align_equ.mgrid');
mgrid_fp = cellstr(mgrid_fp);

for f = 1:size(mgrid_fp,1)
    mgrid2els(mgrid_fp{f});
end

