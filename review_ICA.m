function review_ICA( ICA_results, ICA_infos, SamplingFreq, Segment )
% review_ICA: play with ICA components to look at the effect of their
% removal
%
% review_ICA( ICA_results, Segment )
%
%  Inputs
% --------
% ICA_results : path to ICA results file (.mat), should contain (at least)
% the following variables:
%   - Winv
%   - Activations
% ICA_infos : path to ICA supplementary information files
%
% SamplingFreq: EEG sampling frequency
%
% Segment : indices of time frames to consider (it is not recommended to
% look at the whole reconstruction when it is too long because display and
% refresh might be slow!)
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, November 2018, updated May 2019
%-------------------------------------------------------------------------

% load results and infos
load(ICA_results);
load(ICA_infos);

if nargin > 3
    Activations = Activations(:,Segment); %#ok<NODEF>
end
ThreshFreq = 30;
if SamplingFreq > ThreshFreq
    activationsNEW = zeros(size(Activations,1),length(1:(SamplingFreq/ThreshFreq):size(Activations,2)));
    for c = 1:size(activationsNEW,1)
        % quick and dirty downsampling
        activationsNEW(c,:) = interp1(1:size(Activations,2),Activations(c,:),1:(SamplingFreq/ThreshFreq):size(Activations,2));
    end
    ThisSamplingFreq = ThreshFreq;
    Activations = activationsNEW;
end
if size(Activations,2) > 20000
    warning('Data is big and might take a while to display...')
end

%--------------------------------------------------------------------------
% show EEG with interactive IC removal
%--------------------------------------------------------------------------
plot_ICA_recon(Activations,Winv,ThisSamplingFreq);
title('EEG reconstruction')

%--------------------------------------------------------------------------
% show IC time courses
%--------------------------------------------------------------------------
try
    system([ICAinfos.ActivationPath,' /secondary /minimized &']);
catch
    warning('Could not open SEF file with IC activations, is Cartool installed on system and set as default program to open .sef files?')
end

%--------------------------------------------------------------------------
% show frequency spectrum of each IC
%--------------------------------------------------------------------------
FreqRes = find(FreqList>=1,1)-find(FreqList>=0,1);

% heatmap only within 0 to 10 Hz range (mostly interesting for cardiac)
figure; imagesc(log(FreqPower(1:find(FreqList>=10,1),CompWithCardiacFreq))'); %#ok<NODEF>
set(gca,'xtick',1:FreqRes:find(FreqList>=10,1));
set(gca,'xticklabel',cellstr(num2str((0:10)')));
set(gca,'ytick',1:length(CompWithCardiacFreq));
set(gca,'yticklabel',cellstr(num2str(CompWithCardiacFreq')));
title({'Power (\muV^2/Hz) within 0-10 Hz','of potentially cardiac-related ICs'})
xlabel('Frequency (Hz)');
ylabel('IC #');
colorbar;

% plot also all spectra:
figure;semilogy(FreqList,FreqPower(:,1:size(Winv,2)));
ylabel('Power (\muV^2/Hz)'); xlabel('Frequency (Hz)');
Hf = get(gca,'Children'); grid on;
title('Frequency spectrum of independent components')
clickText(Hf,flip(cellstr(num2str((1:size(Activations,1))'))));

%--------------------------------------------------------------------------
% show ICs topography
%--------------------------------------------------------------------------
try
    system([ICAinfos.TopoPath,' /secondary /minimized &']);
catch
    warning('Could not open ICs topography file, is Cartool installed on system and set as default program to open .sef files?')
end

%--------------------------------------------------------------------------
% show ICs together with bipolar montages of fronto-temporal channels (to
% get EOG-like & EMG-like signals) as well as ECG signal if existing
%--------------------------------------------------------------------------
% EOG-like
try
    system([ICAinfos.EOGlikePath,' /secondary /minimized &']);
catch
    warning('Could not open ICs topography file, is Cartool installed on system and set as default program to open .sef files?')
end
% EMG-like
try
    system([ICAinfos.EMGlikePath,' /secondary /minimized &']);
catch
    warning('Could not open ICs topography file, is Cartool installed on system and set as default program to open .sef files?')
end
% ECG
if isfield(ICAinfos,'LPFabsECGpath')
    try
        system([ICAinfos.LPFabsECGpath,' /secondary /minimized &']);
    catch
        warning('Could not open ICs topography file, is Cartool installed on system and set as default program to open .sef files?')
    end
end

%--------------------------------------------------------------------------
% Show "Z"-scores of FASTER metrics on ICs
%--------------------------------------------------------------------------
figure('name','FASTER metrics "Z"-scores');
subplot(1,3,1);
plot(ICAinfos.FASTER.Zscores(:,1));
hold on; plot(find(ICAinfos.FASTER.Zscores(:,1)>3),ICAinfos.FASTER.Zscores(ICAinfos.FASTER.Zscores(:,1)>3,1),'ro');
title('Mean gradient value');
subplot(1,3,2);
plot(ICAinfos.FASTER.Zscores(:,2));
hold on; plot(find(ICAinfos.FASTER.Zscores(:,2)>3),ICAinfos.FASTER.Zscores(ICAinfos.FASTER.Zscores(:,2)>3,2),'ro');
title('Spatial kurtosis');
subplot(1,3,3);
plot(ICAinfos.FASTER.Zscores(:,3));
hold on; plot(find(ICAinfos.FASTER.Zscores(:,3)>3),ICAinfos.FASTER.Zscores(ICAinfos.FASTER.Zscores(:,3)>3,3),'ro');
title('Hurst exponent');
% legend('Mean gradient value','Spatial kurtosis','Hurst exponent');
% title('FASTER metrics "Z"-scores')

%--------------------------------------------------------------------------
% show correlation of components with EOG / EMG / ECG
%--------------------------------------------------------------------------
if isfield(ICAinfos.ICAcorr,'CompECGcorr')
    figure;
    subplot(1,3,1); imagesc(ICAinfos.ICAcorr.CompEOGcorr); colorbar; title({'Components and bipolar EEG','montage for ocular artefacts'});
    subplot(1,3,2); imagesc(ICAinfos.ICAcorr.CompEMGcorr); colorbar; title({'Components and bipolar EEG','montage for muscular artefacts'});
    subplot(1,3,3); imagesc(ICAinfos.ICAcorr.CompECGcorr); colorbar; title('Components and ECG signal');
else
    figure;
    subplot(1,2,1); imagesc(ICAinfos.ICAcorr.CompEOGcorr); colorbar; title({'Components and bipolar EEG','montage for ocular artefacts'});
    subplot(1,2,2); imagesc(ICAinfos.ICAcorr.CompEMGcorr); colorbar; title({'Components and bipolar EEG','montage for muscular artefacts'});
end

end
