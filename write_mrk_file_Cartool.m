function write_mrk_file_Cartool(FileName,Time1,Time2,MarkerName)
% write_mrk_file_Cartool writes markers to .mrk file for Cartool
%
% USAGE
% write_mrk_file_Cartool(FileName,T1new,T2new,Namenew)
% 
% INPUTS
% FileName: string, filename WITH extension (.mrk)
% Time1: integer, time beginning
% Time2: integer, time end
% MarkerName: cell array of strings with marker names
% 
% Renaud Marquis @FunctionalBrainMapping, 20.06.2017
% last update 2019-05-15

Time1 = Time1(:);
Time2 = Time2(:);
MarkerName = MarkerName(:);

CartoolMagicNumber = 12;
fileID = fopen(FileName,'w');
fprintf(fileID,'%3s\r\n','TL02');
for n = 1:size(Time1,1)
    formatSpec = [repmat(' ',1,CartoolMagicNumber-size(num2str(Time1(n)),2)),'%d\t',repmat(' ',1,CartoolMagicNumber-size(num2str(Time2(n)),2)),'%d\t%s\r\n'];
    fprintf(fileID,formatSpec,Time1(n),Time2(n),MarkerName{n});
end
fclose(fileID);

end
