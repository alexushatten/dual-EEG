function [hdr,data] = fbmlab_read_is(openfilename)
% load_is: opens a Cartool results of inverse solution data file (.is)
%
% inputs: full path and name of the file to open
%
% outputs:
%   hdr: structure with the following fields
%     .magic4 (char array 'IS03'): a tag from Cartool
%     .numelectrodes (int32): number of electrodes
%     .numsolutionpoints (int32): number of solution (source) points
%     .numregularizations (int32): number of regulartizations (Inv. matrices)
%     .isinversescalar (logical): 
%               1 if scalar values (1D)
%               0 if vector values (3D, xyz)
%   data: 2-D float matrix of values...
%     ...if inverse is scalar, numsolutionpoints-by-numelectrodes
%     ...if inverse is vectorial, numtimeframes-by-3*numsolutionpoints
%     where the x, y and z values of the first solution point are listed
%     first (x1 y1 z1), followed by x2, y2, z2 etc.
%     Matrix = numrows x numcols = numsolutionpoints x numelectrodes
%
% Cartool: http://brainmapping.unige.ch/Cartool.htm
%
% author of this script: martin.seeber@unige.ch
%
%% ATTENTION - THIS is currently a 'beta release' type of function
%% FBM Lab internal use only, report if you find bugs
%
% Renaud Marquis @ FBMlab: "data" is rather numsolutionpoints x numelectrodes x numregularizations if scalar
% and 3*numsolutionpoints x numelectrodes x numregularizations if vectorial
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------


fid=fopen(openfilename);

% read in header
hdr.magic4=char(fread(fid,4,'char')');
hdr.numelectrodes=fread(fid,1,'int32');
hdr.numsolutionpoints=fread(fid,1,'int32');
hdr.numregularizations=fread(fid,1,'int32');
hdr.isinversescalar=logical(fread(fid,1,'char'));

for chan = 1:hdr.numelectrodes
    hdr.ElectrodeNames{chan} = char(fread(fid,32,'char')');
end

for sol = 1:hdr.numsolutionpoints
    hdr.SolutionPointNames{sol} = char(fread(fid,16,'char')');
end

for regval = 1:hdr.numregularizations
    hdr.RegularizationValues(regval) = fread(fid,1,'double');
end

for regname = 1:hdr.numregularizations
    hdr.RegularizationNames{regname} = char(fread(fid,32,'char')');
end

clearvars chan sol regval regname

if hdr.isinversescalar
    data = fread(fid,'float');
    data = reshape(data,hdr.numelectrodes,hdr.numsolutionpoints,hdr.numregularizations);
    data = permute(data,[2,1,3]);
else
    data=fread(fid,'float');
    data=reshape(data,hdr.numelectrodes,3*hdr.numsolutionpoints,hdr.numregularizations);
    data = permute(data,[2,1,3]);
end

fclose(fid);
