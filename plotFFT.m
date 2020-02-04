function plotFFT(f,p,fres,flag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % INPUT
%     - f: frequency bins
%     - Power: can be a vector or a horizontal matrix with various power series to be plotted together.
%     - fres: frequ resolution, e.g. 1Hz, will average however fbins
%     - flag: flag for averaging power vectors and calculating standard deviation
% 
% STATUS: final 27.3.2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% refacto, Renaud Marquis @ FBM lab, Sept 2017

% INPUT CHECKS
smooth = 0;
if nargin>2
    if fres>diff(f(1:2))
        smooth=1;
    end
end
av=0;
if nargin > 3
    if ~isempty(flag)
        if size(p,1)==1, error('Cannot average single power series'), end
        av=1;
    end
end

% Smooth
% average frequency bins
if smooth
    f2=[]; p2=[];
    npta=floor(fres/diff(f(1:2)));
    for i=1:npta:length(f)-npta
        f2=[f2 max(f(i:i+npta))];
        p2=[p2 mean(p(:,i:i+npta),2)];
    end
    f=f2;
    p=p2;
end

% PARAMETERS
% Frequency bands
fb = [1 4 8 12 30 90]; % frequ band
nfb = numel(fb);

% Colors
if av
        col= flipud(cbrewer('qual','Paired',2)); %dark blue, then light
else    col=[0 0 0; cbrewer('seq','Reds',size(p,1)-1) ];
end

% Figure
% figure %#RM,2017-09-29: done outside of plotFFT
vm=max(p(:));
% ax = axes('Parent',gcf,'YScale','log','YMinorTick','on','box','on'); %#RM,2017-09-29
set(gca,'YMinorTick','on','box','on'); %#RM,2017-09-29
hold on
% PLOT
if av
    shadedErrorBar(f,p,{@median,@std},'lineprops','b'); 
%     NB: cannot have shadow with logarithmic plot!
%     set(h.mainLine,'color',col(1,:))
%     set(h.patch,'facecolor',col(2,:))
else
%     y=semilogy(f,p'); %#RM,2017-09-29
    y=plot(f,p'); %#RM,2017-09-29
    for i=1:size(p,1)
        set(y(i),'color',col(i,:),'linewidth',1)
    end
end
% set(ax,'XTick',fb) %#RM, 2017-09-25
set(gca,'XTick',fb) %#RM, 2017-09-25

TimesMaxYLim = 1.3;

Ftemp = figure; semilogy(f,exp(p)'); %#RM, 2017-09-25
ylim([prctile(exp(p(:)),10) (TimesMaxYLim+0.5)*exp(vm)]); %#RM, 2017-09-25
NewYTickLabel = get(get(Ftemp,'Children'),'YTickLabel'); %#RM, 2017-09-25
NewYTick = get(get(Ftemp,'Children'),'YTick'); %#RM, 2017-09-25
close(Ftemp); %#RM, 2017-09-25

set(gca,'YTick',log(NewYTick)); %#RM, 2017-09-25
set(gca,'YTickLabel',NewYTickLabel); %#RM, 2017-09-25
% xlim([0 45]) %#RM, 2017-09-25
xlim([0 100]) %#RM, 2017-09-25
% ylim([prctile(p(:),10) 1.5*vm]) %#RM, 2017-09-25
ylim([prctile(p(:),10) TimesMaxYLim*vm]) %#RM, 2017-09-25
set(gca,'YMinorTick','off') %#RM, 2017-09-25
hold on 

% Greying
for i=1:2:nfb-1 
    h = fill([repmat(fb(i),1,2) repmat(fb(i+1),1,2)], [ylim fliplr(ylim)],[0.85 0.85 0.85],'edgecolor','none'); 
    uistack(h,'bottom')
end

% COSMETICS
xlabel('Frequency (Hz)','fontsize',16,'fontweight','bold')
ylabel('Power (µV^2/frequ)','fontsize',16,'fontweight','bold')
flabels = {'\delta';'\theta';'\alpha';'\beta';'\gamma'};
for i=1:nfb-1
    text(mean([fb(i) fb(i+1)]),1.2*vm,flabels{i},'fontsize',16,'fontweight','bold', 'horizontalalignment','center')
end

%%% SAVE
% saveas(gcf,'Power_spectrum_W_S.pdf')