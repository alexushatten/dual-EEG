function [EEG, Time] = brainstorm_resample(EEG, Time, NewFreq, Method)
% brainstorm_resample: resample (e.g. downsample / decimate) signal to
% desired frequency
%
% [TF, time_out] = process_resample(EEG, Time, NewFreq, Method)
%
%  Inputs
% --------
%     - EEG     : Signal to process [nChannels x nTime]
%     - Time    : [1 x nTime], time values in seconds
%     - NewFreq : New sampling frequency (Hz)
%     - Method = 'resample'
%                'resample-rational'
%                'resample-cascade'
%                'interp-decimate-cascade'
%                'fft-spline'
%
%  Outputs
% ---------
%     - EEG   : Resampled signal
%     - Time : New time vector
%
%-------------------------------------------------------------------------
% NB: if using cascaded-integrator comb (CIC) filter using 'resample-cascade',
% compensation filter will be applied by resample.m automatically. However,
% low-pass filter should be applied before using this function in order to
% avoir aliasing. The cutoff frequency of the low-pass filter (LPF) should be
% NewFreq/2 (Nyquist frequency of targeted sampling rate.
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, May 2019
%-------------------------------------------------------------------------

% Default method
if (nargin < 4) || isempty(Method)
    if exist('fir2', 'file'); % Check if Signal Processing Toolbox is installed
        Method = 'resample-cascade';
    else
        Method = 'fft-spline';
    end
end
% Check output frequency
OldFreq = 1 ./ (Time(2) - Time(1));
if (abs(NewFreq - OldFreq) < 0.05)
    EEG = [];
    return;
end
% Round old frequency at x100
OldFreq = round(OldFreq * 100) / 100;


% ===== RESAMPLE DATA =====
% Filtering using the selected method
switch (Method)
    % Bad
    case 'fft-spline'
        % Anti-alias filter
        if (NewFreq < OldFreq)
            % x = process_bandpass('Compute', x, 256, [], 128 * NewFreq / OldFreq, 'bst-fft-fir', 1);
            EEG = process_bandpass('Compute', EEG, 256, [], 128 * NewFreq / OldFreq);  % Replaced by FT, 16-Jan-2016
        end
        % Spline interpolation
        nbnewpoints  = size(EEG,2) * NewFreq / OldFreq;
        nbnewpoints2 = ceil(nbnewpoints);
        lastpointval = size(EEG,2) / nbnewpoints * nbnewpoints2;
        XX = linspace( 1, lastpointval, nbnewpoints2);
        cs = spline( 1:size(EEG,2), EEG);
        EEG = ppval(cs, XX);
        % SigProc Toolbox: Good but can be slow
    case 'resample'
        % Resample: Signal processing toolbox 'resample'
        EEG = resample(EEG', NewFreq, OldFreq)';
        % SigProc Toolbox: Good but output frequency can be different from what is required
    case 'resample-rational'
        % Resample parameters
        [P,Q] = rat(NewFreq / OldFreq, .0001);
        % Resample
        EEG = resample(EEG', P, Q)';
        % Compute output sampling frequency
        NewFreq = P / Q * OldFreq;
        % SigProc Toolbox: Good and fast
    case 'resample-cascade'
        % Resample: Signal processing toolbox 'resample' (cascade)
        EEG = ResampleCascade(EEG, NewFreq, OldFreq, 'resample');
        % SigProc Toolbox: Not so good, not so fast
    case 'interp-decimate-cascade'
        % Resample: Signal processing toolbox 'interp' + 'decimate' (cascade)
        EEG = ResampleCascade(EEG, NewFreq, OldFreq, 'decimate');
end

% Compute new Time vector
% Time = linspace(Time(1), Time(end), size(x,2));
Time = Time(1) + (0:(size(EEG,2)-1)) ./ NewFreq;
end


%% ========================================================================================
%  ====== RESAMPLING FUNCTIONS ============================================================
%  ========================================================================================

%% ====== RESAMPLE-CASCADE =====
% USAGE: [x,Pfac,Qfac] = process_resample('ResampleCascade', x, NewRate, OldRate, Method)
% INPUT:
%     - x       : Signal to process [nChannels x nTime]
%     - NewRate : New sampling frequency (Hz)
%     - OldRate : Original sampling frequency (Hz)
%     - Method  : 'resample' or 'decimate'
% OUTPUT:
%     - x    : Resampled signal
%     - Pfac : Array of successive upsampling factors
%     - Qfac : List of successive downsampling factors
% NOTE: Requires Signal Processing Toolbox
% AUTHOR: John Mosher, 2010
function [x,Pfac,Qfac] = ResampleCascade(x,NewRate,OldRate,Method)
% Default method: 'resample'
if (nargin < 4)
    Method = 'resample';
end
% Common factors
[P,Q] = rat(NewRate/OldRate);
% We want to upsample by P and downsample by Q to achieve the new rate
% But big numbers cause problems.
Pfac = factor(P);
Qfac = factor(Q);
% Longest number of factors
iFacs = max(length(Pfac),length(Qfac));
% Pad the shorter one to have unity factors
Pfac((length(Pfac)+1):iFacs) = 1;
Qfac((length(Qfac)+1):iFacs) = 1;

% So now we have two factorization lists of the same length, and
% prod(Pfac) / prod(Qfac) = P/Q.
Pfac = sort(Pfac,'descend'); % upsample largest first
Qfac = sort(Qfac,'ascend'); % downsample smallest rates first
Rates = Pfac./Qfac;  % rates per step
CRate = cumprod(Rates); % cumulative resampling rates

% We can't go below min(1,P/Q) without losing information. Because of low-pass filtering, don't be too precise
Problem = CRate < (0.9 * P/Q);
if any(Problem)
    fprintf(1, 'RESAMPLE> Warning: Desired rate is %.f\n', P/Q);
end
if any(Pfac > 10)
    disp(['RESAMPLE> Warning: Upsampling by more than 10 in the cascades, P = ' sprintf('%d ', Pfac)]);
end
if any(Qfac > 10)
    disp(['RESAMPLE> Warning: Downsampling by more than 10 in the cascades, Q = ' sprintf('%d ', Qfac)]);
end

% ===== RESAMPLING =====
switch Method
    % Decimate/interp inputs cannot be vectorized
    case 'decimate'
        % Initialize output parameters
        len_resmp = ceil(size(x,2) * prod(Pfac) / prod(Qfac));
        nRow = size(x,1);
        x_resmp = zeros(nRow, len_resmp);
        % Loop on factors and rows
        for iRow = 1:size(x,1)
            x_tmp = x(iRow,:);
            for i = 1:iFacs
                x_tmp = decimate(interp(x_tmp, Pfac(i)), Qfac(i));
            end
            x_resmp(iRow,:) = x_tmp;
        end
        x = x_resmp;
        % Resample takes vectorized inputs
    case 'resample'
        for i = 1:iFacs
            x = resample(x', Pfac(i), Qfac(i))';
        end
end
end
