function read_write_mrk_file_Cartool(FileName,Time1,Time2,MarkerName)
% read_write_mrk_file_Cartool writes markers to .mrk file for Cartool but
% does not overwrite existing content (preprends new content instead)
%
% USAGE
% read_write_mrk_file_Cartool(FileName,T1new,T2new,Namenew)
% 
% INPUTS
% FileName: string, filename WITH extension (.mrk)
% Time1: integer, time beginning
% Time2: integer, time end
% MarkerName: cell array of strings with marker names
% 
% Renaud Marquis @FunctionalBrainMapping, 20.06.2017

T1 = [];
T2 = [];
ML = {};
if exist(FileName,'file')
    try
        [T1,T2,ML] = read_mrk_file(FileName);
    catch me %#ok<NASGU>
        [T1,T2,ML] = read_mrk_Cartool(FileName);
    end
end

CartoolMagicNumber = 12;

if ~isempty(T1) || ~isempty(T2) || ~isempty(ML)
    Time1 = [Time1;T1];
    Time2 = [Time2;T2];
    MarkerName = [MarkerName;ML];
end

fileID = fopen([FileName],'w');
fprintf(fileID,'%3s\r\n','TL02');

for n = 1:size(Time1,1)
    formatSpec = [repmat(' ',1,CartoolMagicNumber-size(num2str(Time1(n)),2)),'%d\t',repmat(' ',1,CartoolMagicNumber-size(num2str(Time2(n)),2)),'%d\t%s\r\n'];
    fprintf(fileID,formatSpec,Time1(n),Time2(n),MarkerName{n});
end
fclose(fileID);

end
