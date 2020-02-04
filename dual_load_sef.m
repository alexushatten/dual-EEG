function [data des] = dual_load_sef(filename)
% [data des] = load_sef(filename)
%
% INPUT
%  filename: file to load
%
% OUTPUTS
%  data: EEG data (nbchannels x nbsamples)
%  des: info structures
%
% G. Birot, FBMLab, Feb. 2013
% refacto, R Marquis, 2017-07-10, 17:52

ELEC_MAX_LENGTH = 8;

fid = fopen(filename,'r');

% header
des.version = fread(fid,1,'int32');
des.nbchannels = fread(fid,1,'int32');
des.Naux = fread(fid,1,'int32');
des.nbsamples = fread(fid,1,'int32');
des.samplingfreq = fread(fid,1,'float32');
des.year = fread(fid,1,'short');
des.month = fread(fid,1,'short');
des.day = fread(fid,1,'short');
des.hour = fread(fid,1,'short');
des.minute = fread(fid,1,'short');
des.second = fread(fid,1,'short');
des.msecond = fread(fid,1,'short');

des.channelnames = cell(1,des.nbchannels);

for n = 1:des.nbchannels % electrodes
    c = fread(fid,ELEC_MAX_LENGTH,'char');
    ind = find(c == 0);
    if (~isempty(ind))
        c = c(1:(ind(1)-1));
    end
    des.channelnames{n} = char(c(:)');
end

% data
data = fread(fid,des.nbchannels*des.nbsamples,'float32');
% data = reshape(data,des.nbsamples,des.nbchannels)';
data = reshape(data,des.nbchannels,des.nbsamples);
fclose(fid);

end