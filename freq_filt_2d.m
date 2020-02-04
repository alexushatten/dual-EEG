function [ImLP,ImHP] = freq_filt_2d( Im, CutOffType, CutOffVal, Show, Hist )
%freq_filt_2d: applies spatial frequency filter to image at specified
%cut-off to separate low- from high-frequency components
%
% Usage:
%-------------------------------------------------------------------------
% 
% [ImLP,ImHP] = freq_filt_2d( Im, CutOffType, CutOffVal, Show, Hist )
%
% Inputs:
%-------------------------------------------------------------------------
%
% Im: n x n array of the image
%
% CutOffType: 'percentile' or 'value' (raw value), default = 'percentile'
%
% CutOffValue: percentile or value at which to cut, default = 99.9
%
% Show: whether to show or not input and filtered image in figure, default
%       = true
%
% Hist: plot an histogram of frequency contents of the input image, default
%       = false
%
% Outputs:
%-------------------------------------------------------------------------
%
% ImLP: n x n array of input image low-pass filtered at cut-off value
%
% ImHP: n x n array of input image high-pass filtered at cut-off value
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, March 2018
%-------------------------------------------------------------------------

if nargin < 2
    CutOffType = 'percentile';
else
    if (strcmpi(CutOffType,'percentile')+strcmpi(CutOffType,'value'))~=1,error('Unrecognized cut-off type, please use ''percentile'' or ''value'''),end;
end
if nargin < 3
    CutOffVal = 99.9;
else
    if strcmpi(CutOffType,'percentile') && (CutOffVal>=100 || CutOffVal<=0),error('Please specify percentile within between 0 and 100 %'),end
end

if nargin < 4
    Show = true;
end

if nargin < 5
    Hist = false;
end

J = dct2(Im);
J2 = J;
switch lower(CutOffType)
    case 'percentile'
        J2(abs(J2) < exp(prctile(log(abs(J2(:))),CutOffVal))) = 0;
    case 'value'
        J2(abs(J2) < exp(CutOffVal)) = 0;
end
ImLP = idct2(J2);

if Show
    figure; subplot(3,1,1); imshow(Im,[]); colorbar; title('Original image')
    subplot(3,1,2); imshow(ImLP,[]); colorbar; title(['Low-pass filtered image at ',num2str(CutOffVal),' %'])
end

J2 = J;
switch lower(CutOffType)
    case 'percentile'
        J2(abs(J2) > exp(prctile(log(abs(J2(:))),CutOffVal))) = 0;
case 'value'
        J2(abs(J2) > exp(CutOffVal)) = 0;
end
ImHP = idct2(J2);

if Show
    subplot(3,1,3); imshow(ImHP,[]); colorbar; title(['High-pass filtered image at ',num2str(CutOffVal),' %'])
end

if Hist
    figure; hist(log(abs(J(:))),30);
    title('Discrete cosine transform (log) of original image')
end

end

