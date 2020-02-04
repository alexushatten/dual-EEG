% pop_repchan() - Replace bad channel(s) using spherical spline
%                 interpolation
%
% Usage:
%   >> [EEG, com] = pop_repchan(EEG); % pop-up window mode
%   >> [EEG, com] = pop_repchan(EEG, 'parameter1', value1, ...
%                                    'parameter2', value2, ...
%                                    'parametern', valuen);
%
% Inputs:
%   EEG       - EEGLAB EEG structure
%   'chans'   - scalar or vector channel(s) to replace
%
% Optional inputs:
%   'nTerms'  - scalar int > 0 number of terms {default 50}
%   'm'       - scalar int > 1 m {default 4}
%   'lookup'  - string file with channel location coordinates
%   'nFrames' - scalar int > 0 number of frames to interpolate
%               vectorized (trade-off between speed and memory usage)
%               {default 1000}
%
% Outputs:
%   EEG       - EEGLAB EEG structure
%   com       - history string
%
% Note:
%   Channel coordinates should be located on the surface of a (unit)
%   sphere.
%
% Author: Andreas Widmann, University of Leipzig, 2006
%
% See also:
%   sserp(), sserpgfcn(), sserpweights(), unitsph(), pop_chanedit()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2006 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Id$

% refacto and stripped down, RM@FBMlab, May 2019
% NB: uses spherical splines, not 3D splines

function EEG = interpolate_bad_EEG_channels(EEG, XYZ, ToInterpolate, nFrames)

warning('Using spherical (not 3D) splines')
if nargin < 4
    nFrames = 1000;
end
% com = '';

% Unit sphere
E = unitsph(XYZ);

% Channel location coordinate matrices E and F
F = E(ToInterpolate, :);
E = E(setdiff(1:size(EEG,1),ToInterpolate), :);
chanArray = setdiff(1:size(EEG,1),ToInterpolate);

% G matrix
[Ginv, g] = sserpgfcn(E, F, 'sp', 0);

% % Reshape EEG.data
% if EEG.trials > 1
%     EEG.data = reshape(EEG.data, EEG.nbchan, []);
% end

% Initialize progress bar
nProgBarSteps = 20;
progBarArray = ceil(linspace(size(EEG, 2) / nProgBarSteps, size(EEG, 2), nProgBarSteps));
progBarHandle = waitbar(0, '0% done', 'Name', 'Replace bad channel(s) -- interpolate_bad_EEG_channels()');
tic

blockArray = [1 : nFrames : size(EEG, 2) size(EEG, 2) + 1];
for iBlock = 1 : length(blockArray) - 1

    % C matrix
    C = sserpweights(EEG(chanArray, blockArray(iBlock) : blockArray(iBlock + 1) - 1), Ginv);

    % Interpolation
    EEG(ToInterpolate, blockArray(iBlock) : blockArray(iBlock + 1) - 1) = sserp(C, g, 'sp');

    % Update progress bar
    if blockArray(iBlock + 1) - 1 >= progBarArray(1)
        progBarArray(1) = [];
        p = (nProgBarSteps - length(progBarArray)) / nProgBarSteps;
        waitbar(p, progBarHandle, [num2str(p * 100) '% done, ' num2str(ceil((1 - p) / p * toc)) ' s left']);
    end

end

% Deinitialize progress bar
if exist('progBarHandle', 'var')
    close(progBarHandle)
end

% % Reshape EEG.data
% if EEG.trials > 1
%     EEG.data = reshape(EEG.data, EEG.nbchan, EEG.pnts, EEG.trials);
% end

% % History string
% com = [inputname(1) ' = interpolate_bad_EEG_channels(' inputname(1) ', ' arg2str(Arg) ');'];

end
