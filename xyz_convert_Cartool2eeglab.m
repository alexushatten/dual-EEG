function xyz_convert_Cartool2eeglab( InputFilename, OutputFilename )
% xyz_convert_Cartool2eeglab: convert .xyz file from Cartool to
% EEGLAB-compatible format, by essentially removing the header and adding
% column 1 with electrode number
%
% xyz_convert_Cartool2eeglab( InputFilename, OutputFilename )
%
%  Inputs
% --------
% InputFilename: string, filepath to Cartool .xyz file
% OutputFilename: string, filepath to new .xyz file compatible with eeglab
%   (if not provided, the suffix "_4eeglab" will be appended to the original
%   file)
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, Date
%-------------------------------------------------------------------------

if nargin < 2
    OutputFilename = spm_file(InputFilename,'suffix','_4eeglab');
end

[x,y,z,name] = read_xyz_file(InputFilename);
write_xyz_eeglab(OutputFilename,x,y,z,name);

end

