function write_xyz_eeglab(OutputFilename,x,y,z,name)
% write_xyz_eeglab: write .xyz file compatible with EEGLAB format
%
% Usage:
%-------------------------------------------------------------------------
% 
% write_xyz_eeglab(Filename,x,y,z,name)
%
% Inputs:
%-------------------------------------------------------------------------
%
% OutputFilename: filename for output file
%
% x, y, z coordinates
%
% name: name of electrode / contact
%
% Xyz: the whole text file as a data array
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, August 2018
%-------------------------------------------------------------------------

fileID = fopen(OutputFilename,'w');

if ~all([numel(x),numel(y),numel(z),numel(name)]==numel(x))
    error('Input size mismatch')
end
for l = 1:numel(x)
    formatSpec = '\t%d %.7f %.7f %.7f %s\r\n';
    fprintf(fileID,formatSpec,l,x(l),y(l),z(l),name{l});
end

fclose(fileID);

end
