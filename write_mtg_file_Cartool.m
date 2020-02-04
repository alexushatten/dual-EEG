function write_mtg_file_Cartool(FileName,labelsB2)
% write_mtg_file_Cartool writes montage (.mtg) file for Cartool based on
% electrodes pairs
%
% USAGE
% write_mtg_file_Cartool(FileName,T1new,T2new,Namenew)
% 
% INPUTS
% FileName: string, filename WITH extension (.mtg)
% labelsB2: n x 2 cell array of strings with electrode names
% 
% Renaud Marquis @FunctionalBrainMapping, February 2018

fileID = fopen([FileName],'w');
fprintf(fileID,'%3s\r\n','MT01');
for n = 1:size(labelsB2,1)
    fprintf(fileID,'%s\t%s\r\n',labelsB2{n,1},labelsB2{n,2});
end
fclose(fileID);

end
