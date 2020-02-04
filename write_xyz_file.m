function write_xyz_file(filename,x,y,z,name)
%write_xyz_file: writes Cartool .xyz file
%
% Usage:
%-------------------------------------------------------------------------
% 
% write_xyz_file(filename, x,y,z,name)
%
% Inputs:
%-------------------------------------------------------------------------
%
% x, y, z coordinates
%
% name: name of electrode / contact
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, October 2018
%-------------------------------------------------------------------------

if ~all(~([size(x,1),size(y,1),size(z,1),size(name,1)]-size(x,1)))
    error('Input size mismatch')
end

fileID = fopen(filename,'w');
r = 125.649;
fprintf(fileID,'%d\t%.3f\r\n',size(x,1),r); % the radius is anyhow not used anymore in Cartool so you can put any value for "r"

CartoolMagicNumber = 16;

for l = 1:size(x,1)
    FormatSpec = ['  %.7e  %.7e  %7e',repmat(' ',1,CartoolMagicNumber-size(name{l},2)),'%s\r\n'];
    fprintf(fileID,FormatSpec,x(l),y(l),z(l),name{l});
end

fclose(fileID);

end
