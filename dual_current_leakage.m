function [Pch,Pngb,f,ngb,Ratio,psigntest,psignrank,pttest]=dual_current_leakage(EEG,fs,ic_hd_montage,Method)

% INPUTS
% - EEG: hdEEG with 257 channels 
% - fs: frequency sampling (sampling rate)
% - ic_hd_montage: correspondance of scalp and stereo electrodes
%   - 1st column: name of stereo electrode
%   - 2nd column: neighbours of stereo electrode
%   - 3rd column: closest one, two or three electrodes, entered manually
%   based on photos
% - Method: 3D distance with 25 mm as threshold ('3d') or neighbours part
% of the hexagone around the electrode ('hex'). In the latter case,
% electrodes #61 and #220 have 7 and not six neighbours, and the threshold
% distance is 26.75 mm instead.
%
% Maxime Beaud
%
% refacto, Renaud Marquis @ FBM lab

% hdEEG channels close to stereo
chs = [];
for i=1:size(ic_hd_montage,1)
    chs=[chs; ic_hd_montage{i,3}(:)];
end
chs(isnan(chs))=[]; %if stereo electrode not close to any hdEEG ch

switch Method
    case '3d'
        % hdEEG channels around those channels
        chi = [];
        ngb=[]; % neighbours
        for i=1:numel(chs)
            newngb = neighbours_256chs(chs(i));
            newngb = setdiff(newngb,chs); % do not take chs
            ngb = [ngb; newngb(:)];
            chi = [chi ; repmat(chs(i),numel(newngb),1)];
        end
    case 'hex'
        % ALTERNATIVE (#RM,2017-09-29):
        chi = [];
        ngb=[]; % neighbours
        for i=1:numel(chs)
            newngb = my_neighbours_256chs(chs(i))';
            newngb = setdiff(newngb,chs); % do not take chs
            ngb = [ngb; newngb(:)];
            chi = [chi ; repmat(chs(i),numel(newngb),1)];
        end
end

ngb = [ngb, chi];

% POWER in chs
% Pch=nan(numel(chs),501);
Pch=nan(numel(chs),2*fs+1);
for i = 1:numel(chs)
    [p,f]=FFTeeg(EEG(chs(i),:),fs,4,.5,0);
    Pch(i,:)=mean(p);
end


% POWER in neighbours
% Pngb=nan(size(ngb,1),501);
Pngb=nan(size(ngb,1),2*fs+1);
for i = 1:size(ngb,1)
    [p]=FFTeeg(EEG(ngb(i,1),:),fs,4,.5,0);
    Pngb(i,:)=mean(p);
end

Pch = log(Pch); %#RM,2017-09-29
Pngb = log(Pngb); %#RM,2017-09-29

% RATIO
% Frequency bands
fb = [1 4 8 12 30 90]; % frequ band
nfb = numel(fb)-1;
r = nan(numel(chs),nfb);
for i = 1:nfb
    for j=1:numel(chs)
        fidx = f>=fb(i) & f<= fb(i+1);
        pch = mean(Pch(j,fidx));
        pngb = mean(mean(Pngb(ngb(:,2)==chs(j),fidx)));
        r(j,i) = pch/pngb; 
    end
end
Ratio = r;
[~,pttest]=ttest(r,1); % t-test on mean
psigntest = nan(1,(numel(fb)-1));
for i = 1:(numel(fb)-1)
    psigntest(i)=signtest(r(:,i),1); % Wilcoxon rank sum test (because even log of power spectral density is not at all normally distributed)
end
psignrank = nan(1,(numel(fb)-1));
for i = 1:(numel(fb)-1)
    psignrank(i)=signrank(r(:,i),1); % Wilcoxon signed rank test (are different but neighbouring electrodes in the same subject sampled at the same moment different samples (entities) or not? (they share at least electrical sources... but can we really consider these as repeated measures?)
end
r_sd = std(r);
r = mean(r);

% FIGURE
figure; subplot(2,1,1); %#RM,2017-09-29
fres = 0.5;
plotFFT(f,Pch,fres,1)
hold on
f2=[]; p2=[];
npta=floor(fres/diff(f(1:2)));
for i=1:npta:length(f)-npta
    f2=[f2 max(f(i:i+npta))];
    p2=[p2 mean(Pngb(:,i:i+npta),2)];
end
f=f2;
Pngb=p2;
shadedErrorBar(f,Pngb,{@median,@std},'lineprops','r'); 

% figure %#RM,2017-09-29
subplot(2,1,2); %#RM,2017-09-29
x = 1:nfb;
vm = max(r(:))+max(r_sd(:));
bar(x,r)
hold on
errorbar(x,r,r_sd,'linestyle','none')
xlim([.5 nfb+.5])
ylim([0 1.3*vm])
line(xlim,[1 1],'linestyle','-.') %bsl

xlabel('Frequency band [Hz]','fontsize',16,'fontweight','bold')
ylabel('Relative Power [%]','fontsize',16,'fontweight','bold')
flabels = {'\delta';'\theta';'\alpha';'\beta';'\gamma'};
% set(gca,'xtick',x,'xticklabel',[],'ytick',[.5 1 1.5])
for i=1:nfb
    text(i,1.1*vm,sprintf('p = %1.3f',psignrank(i)),'fontsize',12,'fontweight','bold','horizontalalignment','center')
    text(i,1.2*vm,flabels{i},'fontsize',16,'fontweight','bold','horizontalalignment','center')
end
