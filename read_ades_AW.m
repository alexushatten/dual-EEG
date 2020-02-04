function [Channels, Modality, Nsamples, FreqSamp] = read_ades_AW(Filename)
% read_ades_AW: read AnyWave .ades header file
%
% [Channels, Nsamples, FreqSamp] = read_ades_AW(Filename)
%
%  Inputs
% --------
% Filename: path to AnyWave .ades file
%
%  Outputs
% ---------
% Channels: cell array with channel names
% Modality: cell array with modality for each channel
% Nsamples: number of samples
% FreqSamp: sampling rate
%
%-------------------------------------------------------------------------
% http://meg.univ-amu.fr/wiki/AnyWave:ADES
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, March 2019
%-------------------------------------------------------------------------

fid = fopen(Filename,'r');
Header = fgetl(fid);
if ~strcmp(Header,'#ADES header file ')
    warning('Unexpected header of file %s, please make sure the file is .ades AnyWave metadata file')
end
Count = 1;
while (feof(fid)==0)    % while  not at the end of the file
    line = fgetl(fid); %take line by line
    pieces = regexp(line,'=','split');
    if strcmp(pieces{1},'samplingRate ')
        FreqSamp = str2double(pieces{2});
    elseif strcmp(pieces{1},'numberOfSamples ')
        Nsamples = str2double(pieces{2});
    else
        if ~isempty(pieces{1})
            Channels{Count} = strtrim(pieces{1}); %#ok<AGROW>
            Modality{Count} = strtrim(pieces{2}); %#ok<AGROW>
            Count = Count+1;
        end
    end
    
end
fclose(fid);
Channels = Channels(:);
Modality = Modality(:);

end
