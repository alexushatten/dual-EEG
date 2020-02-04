function save_is(filename, is)
% save_is: writes a Cartool results of inverse solution data file (.is)
%
%   filename: full path of file to write
%
%   is: structure with following fields:
%       .magic4 (char array 'IS03'): a tag from Cartool
%       .numelectrodes (int32): number of electrodes
%       .numsolutionpoints (int32): number of solution (source) points
%       .numregularizations (int32): number of regulartizations (Inv. matrices)
%       .isinversescalar (logical):
%             1 if scalar values (1D)
%             0 if vector values (3D, xyz)
%       .ElectrodeNames : electrode names
%           (cell of strings, 1 x is.numelectrodes with each name being a
%               1 x 32 characters-long string)
%       .SolutionPointNames : solution point names
%           (cell of strings, 1 x is.numsolutionpoints with each name being
%               a 1 x 16 characters-long string)
%       .RegularizationValues : regularization values (double, usually 1 x 13)
%       .RegularizationNames : regularization names
%           (cell, usually 1 x 13 with "Regularization XX" where XX is the
%               column number, with name being a 1 x 32 characters-long string)
%       .data: 2-D float matrix of values...
%         ...if inverse is scalar, numtimeframes-by-numsolutionpoints
%         ...if inverse is vectorial, numtimeframes-by-3*numsolutionpoints
%         where the x, y and z values of the first solution point are listed
%         first (x1 y1 z1), followed by x2, y2, z2 etc.
%         Matrix = numrows x numcols = numsolutionpoints x numelectrodes
%
% WARNING /!\ currently only 'IS03' inverse solution matrix types are
% supported!! 'ISO1' and 'ISO2' are not commonly used anymore
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, May 2018, based on scripts from martin.seeber@unige.ch & G. Birot
%-------------------------------------------------------------------------

fid=fopen(filename, 'w');

% header
fwrite(fid, 'IS03', 'char');
fwrite(fid, is.numelectrodes, 'int32');
fwrite(fid, is.numsolutionpoints, 'int32');
fwrite(fid, is.numregularizations, 'int32');
fwrite(fid, is.isinversescalar, 'char');

for chan = 1:is.numelectrodes
    % force 32-characters-long string
    fwrite(fid, pad(is.ElectrodeNames{chan}(:)',32),'char');
end

for sol = 1:is.numsolutionpoints
    % force 16-characters-long string
    fwrite(fid, pad(is.SolutionPointNames{sol}(:)',16),'char');
end

for regval = 1:is.numregularizations
    fwrite(fid, is.RegularizationValues(regval), 'double');
end

for regname = 1:is.numregularizations
    % force 32-characters-long string
    fwrite(fid, pad(is.RegularizationNames{regname}(:)',32),'char');
end

% data
Data = permute(is.data,[2 1 3]);
fwrite(fid, Data(:), 'float32');

fclose(fid);

end
