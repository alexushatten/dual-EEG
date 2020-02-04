function h = plot_spike( EEG, ChanMax, Srate, ChanelLabels, Offset, EEGoverlay )
% plot_spike: plot spike
%
% plot_spike( EEG, Srate, ChanelLabels, Offset )
%
%  Inputs
% --------
% EEG : t x c EEG traces (matrix where t is time frame and c is channel)
%
% ChanMax : channel to highlight in blue
%
% Srate : sampling rate in Hz
%
% ChanelLabels: 1 x c cell array of strings with channel labels
%
% Offset: offset between neighbouring channel traces, default is 2x range
% from 5th to 95th percentile [5 95]
%
% EEGoverlay: other EEG trace to overlay on top of the primary EEG trace
% (optional, same format)
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, July 2018, last change April 2019
%-------------------------------------------------------------------------

ChanelLabels = ChanelLabels(:);

if nargin < 5
    Offset = [5 95];
end

h = figure('units','normalized','outerposition',[0.1 0.1 0.4 0.75],'Visible','off');
hold on;

PixSS = get(0,'ScreenSize'); PixSS(3) = PixSS(3)/1.8750;
if nargin > 5
    EEGrange = prctile([EEG,EEGoverlay],Offset(2),1)-prctile([EEG,EEGoverlay],Offset(1),1);
    ofst = (1:size(EEG,2))*nanmean(EEGrange)*2 + 0.001;
    EEGp = bsxfun(@plus, EEG, fliplr(ofst))';
    EEGoverlayp = bsxfun(@plus, EEGoverlay, fliplr(ofst))';
    plot(EEGp','k');
    if ~isempty(ChanMax)
        plot(EEGp(ChanMax,:),'b');
    end
    plot(EEGoverlayp','b');
else
    EEGrange = prctile(EEG,Offset(2),1)-prctile(EEG,Offset(1),1);
    ofst = (1:size(EEG,2))*nanmean(EEGrange)*2 + 0.001;
    EEGp = bsxfun(@plus, EEG, fliplr(ofst))';
    plot(EEGp','k');
    if ~isempty(ChanMax)
        plot(EEGp(ChanMax,:),'b');
    end
end

Yrange = [nanmean(EEGp(1,:)),nanmean(EEGp(end,:))];
Ylim = [Yrange(2)-(nanmean(EEGrange)*2) Yrange(1)+(nanmean(EEGrange)*2)];
ylim(Ylim);
xlim([1 1*Srate+1]);

ZoomH = zoom(h);
set(ZoomH,'Motion','horizontal');
PanH = pan(h);
set(PanH,'Motion','horizontal');

if nargin > 3
    set(gca,'ytick',linspace(round(Yrange(2)),round(Yrange(1)),numel(ChanelLabels)));
    set(gca,'yticklabel',flip(ChanelLabels));
else
    set(gca,'ytick',[]);
end

H = get(gca,'Children');

if nargin > 5
    if ~isempty(ChanMax)
        highlightChannel(H,repmat(flip([ChanelLabels;{''}]),1,2));
    else
        highlightChannel(H,repmat(flip(ChanelLabels),1,2));
    end
else
    if ~isempty(ChanMax)
        highlightChannel(H,flip([ChanelLabels;{''}]));
    else
        highlightChannel(H,flip(ChanelLabels));
    end
end

%% GUI buttons

ButtonPositionsX = 0:round(PixSS(3)*0.75)/3:round(PixSS(3)*0.75);
ButtonPositionsX(1)=ButtonPositionsX(1)+round(0.037*PixSS(4));
ButtonPositionsX(2)=ButtonPositionsX(2)-round(0.037*PixSS(4));
ButtonPositionsX(3)=ButtonPositionsX(3)-round(0.037*PixSS(4));
ButtonPositionsX(4)=ButtonPositionsX(4)-round(0.037*2*PixSS(4));

Buttons(1) = uicontrol('Style', 'pushbutton', 'String', 'OK, continue',...
    'Position', [round(mean([ButtonPositionsX(2) ButtonPositionsX(3)])) round(0.0185*PixSS(4)) round(1/9*PixSS(4)) round(0.0185*PixSS(4))],...
    'Callback', @close_EEG);

set(Buttons,'Units','normalized');

function close_EEG(source, callbackdata) %#ok<INUSD>
    h = gcf;
    close(h);
end

h.Visible = 'on';

end

function highlightChannel(h,txt,interpreter)
%function highlightChannel(h,txt,interpreter)
%
% Adds code to a figure object so that when you click on it, a text box
% appears with the desired text in it. When you click on the text box, it
% disappears.
%
% Required Inputs:
%  h   - figure object handle (vector or singleton)
%  txt - string of cell array of strings containing object labels
%
% Optional Inputs:
%  interpreter - The value of the text object's "interpreter" property.
%                {default: 'tex'}
%
% Author: David Groppe
% Mehtalab, 2012
% Renaud Marquis, FBMlab: added thickening of linewidth when clicked
% refacto from clickText

if nargin<3,
    interpreter='tex';
end

hTparams=['set(ht,''backgroundcolor'',''w'',''horizontalalignment'',''center'',''verticalalignment'',''middle'',''interpreter'',''' interpreter ''',''buttondownfcn'',''delete(gcbo);'');'];

if iscell(txt)
    if length(h)~=length(txt)
        error('To the number of elements of h and txt are different.');
    end
    
    for a=1:length(h),
        set(h(a),'userdata',txt{a});
        set(h(a),'linewidth',0.5); % #RM@FBMlab,December 2017
%         set(h(a),'color','k');
        %         bdfcn=['Cp = get(gca,''CurrentPoint''); '  'Xp=Cp(1,1);',  'Yp=Cp(1,2);',   'Zp=Cp(1,3);',  'dat=get(gcbo,''userdata'');',   'ht=text(Xp,Yp,Zp,sprintf(''%s'',dat));',  hTparams]; % #RM@FBMlab,December 2017
        bdfcn=['Cp = get(gca,''CurrentPoint''); '  'Xp=Cp(1,1);',  'Yp=Cp(1,2);',   'Zp=Cp(1,3);',  'dat=get(gcbo,''userdata'');',   'Xlim = xlim; ht=text(round(Xlim(2)/2),-20,-5,sprintf(''%s'',dat)); if get(gco,''linewidth'')==0.5,set(gco,''linewidth'',3),set(gco,''color'',''red''),else set(gco,''linewidth'',0.5),set(gco,''color'',''k''),end; ',  hTparams]; % #RM@FBMlab,December 2017
        set(h(a),'buttondownfcn',bdfcn);
    end
else
    set(h,'userdata',txt);
    set(h,'linewidth',0.5); % #RM@FBMlab,December 2017
%     set(h,'color','k');
    %     bdfcn=['Cp = get(gca,''CurrentPoint''); '  'Xp=Cp(1,1);',  'Yp=Cp(1,2);',   'Zp=Cp(1,3);', 'dat=get(gcbo,''userdata'');', 'ht=text(Xp,Yp,Zp,sprintf(''%s'',dat));',  hTparams]; % #RM@FBMlab,December 2017
    bdfcn=['Cp = get(gca,''CurrentPoint''); '  'Xp=Cp(1,1);',  'Yp=Cp(1,2);',   'Zp=Cp(1,3);', 'dat=get(gcbo,''userdata'');', 'Xlim = xlim; ht=text(round(Xlim(2)/2),-20,-5,sprintf(''%s'',dat)); if get(gco,''linewidth'')==0.5,set(gco,''linewidth'',3),set(gco,''color'',''red''),else set(gco,''linewidth'',0.5),set(gco,''color'',''k''),end; ',  hTparams]; % #RM@FBMlab,December 2017
    set(h,'buttondownfcn',bdfcn);
end

end
