function plot_ERP( EEG, HeightTresh, TrackNames, SigTracks )
% plot_ERP: interactive plot of evoked potential
%
% plot_ERP( EEG, HeightTresh, TrackNames, SigTracks )
%
%  Inputs
% --------
% EEG: [channels x time] evoked potential
% HeightThresh: threshold of significance for shaded rectangle
% TrackNames: cell array of strings with channel labels
% SigTracks: indices of channels with significant response
%
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, November 2018
%-------------------------------------------------------------------------

figure('units','normalized','outerposition',[0.1 0.1 0.75 0.75]);

if nargin > 1
    boundedline(1:1000,zeros(1000,1),repmat(HeightTresh,1,1000),'k'); % outlinebounds(Lb,Ub);
end

% hold on; set(gca,'XTick',1:50:1000); set(gca,'XTickLabel',{-200:50:800});
ylabel('amplitude (uV)'); % ylabel('t-score');
xlabel('time (samples)'); grid on;
plot(EEG');

if nargin > 2
    H = get(gca,'Children');
    clickText(H,flipud(vertcat({'';''},TrackNames(:)))); % the two last are for the zeros (black line) and the boundedline (grey patch)
else
    H = get(gca,'Children');
    clickText(H,flipud(vertcat({'';''},cellstr(num2str([1:size(EEG,1)]'))))); % the two last are for the zeros (black line) and the boundedline (grey patch)
end

if nargin > 3
    NonSigTracks = ~dnif(SigTracks,length(TrackNames));
    set(H(find(flipud(vertcat([0;0],NonSigTracks)))),'linestyle',':'); %#ok<FNDSB>
%     set(H(find(flipud(vertcat([0;0],dnif(SigTracks,length(TrackNames)))))),'linewidth',2); %#ok<FNDSB>
end

end

