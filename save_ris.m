function save_ris(filename, ris)
% function save_ris(filename, ris)
%
% Write RIS Cartool file
%
% inputs
%   filename: file name
%   ris: result of inverser solution structure
%        .nSp: number of solution points
%        .nTimeFrame: number of time frames
%        .samplingFreq: sampling frequency
%        .data: result of inverse [nSp x nTimeFrame] (or [3 x nSp x nTimeFrame])
%
% G. Birot - May 2014
% refacto & added support for vectorial RIS, Renaud Marquis @ FBM lab, August 2018

fid = fopen(filename, 'w');

% header
fwrite(fid, 'RI01', 'char');
fwrite(fid, ris.nSp, 'int32');
fwrite(fid, ris.nTimeFrame, 'int32');
fwrite(fid, ris.samplingFreq, 'float32');
if numel(ris.data)==(ris.nSp*ris.nTimeFrame)
    FlagVect = false;
    fwrite(fid, 1, 'char');
else
    FlagVect = true;
    fwrite(fid, 0, 'char');
end

% data
if ~FlagVect
    fwrite(fid, ris.data, 'float32'); % to verify => #RM: transpose removed and description of required format added in help section
else
%     ris.data = reshape(ris.data,3,ris.nSp,ris.nTimeFrame);
    fwrite(fid, ris.data, 'float32'); % #RM: TO CHECK
end

fclose(fid);

end
