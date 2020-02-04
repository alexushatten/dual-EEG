function [RisMat, IsInverseScalar, Fs] = read_ris_Cartool(Filename)
% read_ris_Cartool: reads Cartool .ris (Results of Inverse Solution computation) file
%
% Usage:
%-------------------------------------------------------------------------
% 
% [RisMat, IsInverseScalar, Fs] = read_ris_Cartool(Filename)
%
% Outputs:
%-------------------------------------------------------------------------
%
% Sources: matrix of sources
%
% IsInverseScalar: whether the Results of Inverse Solution are in scalar
% (norm) or vectorial form
%
% Fs: sampling frequency
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, April 2018, based on functions written by G. Birot (May
% 2014) and Rémi Tyrand, FBM Lab, Hôpitaux Universitaires de Genève (26/10/2011).
%-------------------------------------------------------------------------

fid = fopen(Filename,'r');
Hdr = fread(fid,4,'*char')'; % should be RI01
if ~strcmp(Hdr,'RI01')
    warning('Unexpected header of .ris file (not ''RI01''), please check input file!')
end
Temp = fread(fid,2,'*int32'); Nsp = Temp(1); Ntf = Temp(2);
Fs = fread(fid,1,'*float');
IsInverseScalar = fread(fid,1,'*int8');
RisVec = fread(fid,'*float');
fclose(fid);
if IsInverseScalar
    RisMat = reshape(RisVec,Nsp,Ntf); % counter-intuitive but works
else
    RisMat = reshape(RisVec,Nsp*3,Ntf);
    RisMat = reshape(RisMat',Ntf,3,Nsp); % looks stupid but works
    RisMat = permute(RisMat,[3 1 2]);
end

end
