function Signal = brainstorm_notch(Signal, Fs, FreqList)
% brainstorm_notch: Remove one or more sinusoids from a signal
%
% USAGE:      Signal = brainstorm_notch(Signal, Fs, FreqList)
%
% similar to iirnotch with less dB removed than default -Ab = -3dB

% @=============================================================================
% This function is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2017 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Francois Tadel, 2014-2015
% 
% Code inspired from MatlabCentral post:
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/292960
%
% refacto based on process_notch, Renaud Marquis @ FBM lab, April 2018

% Define a default width
FreqWidth = 1;
% Remove the mean of the data before filtering
xmean = mean(Signal,2);
Signal = bst_bsxfun(@minus, Signal, xmean);
% Remove all the frequencies sequencially
for ifreq = 1:length(FreqList)
    fprintf('Notch %d / %d... \n',ifreq,length(FreqList));
    % Define coefficients of an IIR notch filter
    delta = FreqWidth/2;
    % Pole radius
    r = 1 - (delta * pi / Fs);
    theta = 2 * pi * FreqList(ifreq) / Fs;
    % Gain factor
    B0 = abs(1 - 2*r*cos(theta) + r^2) / (2*abs(1-cos(theta)));
    % Numerator coefficients
    B = B0 * [1, -2*cos(theta), 1];
    % Denominator coefficients
    A = [1, -2*r*cos(theta), r^2];
    % Filter signal
    Signal = filtfilt(B,A,Signal')';
    
end
% Restore the mean of the signal
Signal = bst_bsxfun(@plus, Signal, xmean);

end