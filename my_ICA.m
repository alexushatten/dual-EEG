function [ W, Winv, Activations, Topos, CompWithCardiacFreq, FreqPower, FreqList ] = my_ICA( EEG, SamplingFreq )
% my_ICA: perform ICA decomposition of EEG data. EEG traces are assumed to be
% already filtered, and should not contain bad segments or bad channels.
%
% [ W, Winv, Activations, Topos, ...
%   CompWithCardiacFreq, FreqPower, FreqList ] = my_ICA( EEG, SamplingFreq )
%
%  Inputs
% --------
% EEG: [channels x time] EEG traces
% SamplingFreq: [1 x 1] integer, sampling rate of EEG
%
%  Outputs
% ---------
%                   W : [components x channels] unmixing matrix
%                Winv : [channels x components] mixing matrix
%         Activations : [components x time] time courses of components
%               Topos : [components x channels] topography of components (transposed mixing matrix)
% CompWithCardiacFreq : [1 up to components] indices of components highly correlated with
%                       frequency of cardiac beats
%           FreqPower : [2*SamplingRate+1 x components] frequency power of components time course
%            FreqList : [2*SamplingRate+1 x 1] list of frequency bands corresponding to FreqPower
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, May 2019
%-------------------------------------------------------------------------
% #TODO: call pop_runica.m instead of runica.m, such that it can also perform
%        ICA using other algorithms (binica, JADE, fastica), see for example
%        https://github.com/openroc/eeglab/blob/master/tags/EEGLAB7_0_0_0beta/functions/popfunc/pop_runica.m
%-------------------------------------------------------------------------

tic;
%--------------------------------------------------------------------------
% Calcul Number of Principal Components if data reduction is needed
%--------------------------------------------------------------------------
fprintf('Calculating number of principal components\n');
% k = 20;
% "20" was changed to "50" in one of the last version of FBMlab code, but here we often end up with (way) less than 100 components, which might be not enough!
k = 30; % However, as a general rule, other studies, guidelines and automated pipelines suggest instead to set it to 30, see for example https://doi.org/10.3389/fnins.2018.00097 & https://sccn.ucsd.edu/wiki/Chapter_09:_Decomposing_Data_Using_ICA
NPCA=fix(sqrt(size(EEG,2)/k));
if NPCA<size(EEG,1)
    Nact=NPCA;
    fprintf('Data reduction will be performed before ICA using %d principal components\n\n',NPCA);
else
    Nact=size(EEG,1);
    NPCA=0;
end

%--------------------------------------------------------------------------
% ICA decomposition
%--------------------------------------------------------------------------

fprintf('ICA decomposition\n');
if NPCA
    [weights, sphere] = runica(EEG,'pca',NPCA);
else 
    [weights, sphere] = runica(EEG);    
end

% Calculate unmixing matrix
fprintf('Calculating unmixing matrix\n');
W=weights*sphere;
% clearvars weights sphere;

% Calculate mixing matrix
fprintf('Calculating mixing matrix\n');
Winv=pinv(W);

% Calculate time courses of components
fprintf('Calculating components time-courses\n');
Activations=W*EEG;
% clearvars EEG;

% Calculate topography of components
fprintf('Calculating ICA topographies\n');
Topos=Winv';

%% ----- EXPERIMENTAL -----
%--------------------------------------------------------------------------
% DETECT SUSPICIOUS COMPONENTS
%--------------------------------------------------------------------------
fprintf('Detecting possible ECG-related components\n');
intercorr_max_value=zeros(1,Nact);
% The values below were taken from FBMlab script. It covers a safe range
% of possible heart rates (45 to 150). In the general healthy population,
% heart rate is generally between 60 and 100 at rest, with slightly different
% intervals for childrens.
maxlag=round((60/45)*SamplingFreq);
minlag=round((60/150)*SamplingFreq);
for i=1:Nact
    tmp=xcorr(Activations(i,:),maxlag,'coeff');
    intercorr_max_value(i)=max(tmp(maxlag+minlag:end));
end
CompWithCardiacFreq = find(intercorr_max_value>(mean(intercorr_max_value)+2*std(intercorr_max_value)));
% if ~isempty(cardiac_cpt)
%     cardiac_cpt=cardiac_cpt-1; % Cartool begins at 0
%     fprintf(fid_vrb,'These components may be cardiac-related:\n');
%     for i=1:length(cardiac_cpt)
%         fprintf(fid_vrb,'CP%d, %2.2f%%\n',cardiac_cpt(i),100*intercorr_max_value(cardiac_cpt(i)+1));
%     end
%     fprintf(fid_vrb,'\n');
% end

%--------------------------------------------------------------------------
% frequency analysis of activation vectors
%--------------------------------------------------------------------------
fprintf('Estimating frequency spectrum of components\n');
FreqPower = nan(SamplingFreq*2+1,Nact); % there will be 4 time bins per Hz, so we need 2 times the sampling frequency + 1
% if false % ...DOES NOT WORK RELIABLY BECAUSE OF KASPERSKY ANTIVIRUS AT HUG!! license('test','Distrib_Computing_Toolbox')
%     parfor i = 1:Nact
%         TmpAct = Activations(i,:); % #RM@FBMlab: slices but doesn't broadcast the array (https://stackoverflow.com/questions/46508947/matlab-parfor-loop-freezing-with-a-large-data-set)
% %         tic;[Pxx_p,F_p] = pwelch(TmpAct(:),4*SamplingFreq,0,0:0.5:(round(SamplingFreq/2)),SamplingFreq);toc
%         [p,~] = FFTeeg(TmpAct,SamplingFreq,4,.5,0);toc %figure;semilogy(ff',mean(p,1)); ylabel('Power (\muV^2/Hz)'); xlabel('Frequency (Hz)')
%         % NB: FFTeeg gives results very similar to pwelch() but is about
%         % 5.5 x faster!
%         FreqPower(:,i) = mean(p,1)';
%     end
%     [~,FreqList] = FFTeeg(Activations(1,:),SamplingFreq,4,.5,0); % otherwise FreqList will not be available (because it is internal to the parfor)
% else
    % Fortunately, FFT called by FFTeeg is multithreaded, so the
    % computation will still be relatively fast!
for i = 1:Nact
    TmpAct = Activations(i,:); % #RM@FBMlab: slices but doesn't broadcast the array (https://stackoverflow.com/questions/46508947/matlab-parfor-loop-freezing-with-a-large-data-set)
    [p,FreqList] = FFTeeg(TmpAct,SamplingFreq,4,.5,0); %figure;semilogy(ff',mean(p,1)); ylabel('Power (\muV^2/Hz)'); xlabel('Frequency (Hz)')
    FreqPower(:,i) = mean(p,1)';
end
% end

Toc = toc;
if Toc<60
    fprintf('ICA decomposition (+ frequency of components)\nperformed in %d seconds.\n',round(Toc));
elseif (Toc/60)<60
    fprintf('ICA decomposition (+ frequency of components)\nperformed in %d minutes.\n',round(Toc/60));
else
    fprintf('ICA decomposition (+ frequency of components)\nperformed in %d hours and %d minutes.\n',round(Toc/3600),round(rem(Toc,3600)/60));
end

end

