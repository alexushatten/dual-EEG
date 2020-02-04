function write_lm_file(filename,PathsArray)
%write_xyz_file: writes Cartool .lm (Link Many) file
%
% Usage:
%-------------------------------------------------------------------------
% 
% write_lm_file(filename,PathsArray)
%
% Inputs:
%-------------------------------------------------------------------------
%
% PathsArray: cell of strings containing paths to files to be linked
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, November 2018
%-------------------------------------------------------------------------

[~,~,Ext] = fileparts(filename);
if ~strcmpi(Ext,'.lm')
    error('Please add ".lm" extension to output filename')
end

fileID = fopen(filename,'w');


for l = 1:size(PathsArray,1)
    fprintf(fileID,'%s\r\n',PathsArray{l});
end

fclose(fileID);

end
