function plot_ICA_recon( activations, Winv, Srate, ChanelLabels, Offset )
% plot_ICA_recon: reconstruct EEG with ICA results and play with components
% to interactively look at results when they are removed
%
% plot_ICA_recon( activations, Winv, Srate, ChanelLabels, Offset )
%
%  Inputs
% --------
% 
% activations: [components x time] ICA time courses
%
% Winv: [channels x components] mixing matrix
%
% ChanelLabels: 1 x c cell array of strings with channel labels
%
% Srate: sampling rate of EEG (and ICA time courses)
%
% Offset: offset between neighbouring channel traces, default is 2x range
% from 5th to 95th percentile [5 95]
%
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, November 2018, last updated November 2018
%-------------------------------------------------------------------------

if nargin < 5
    Offset = [5 95];
end
if nargin < 4
    ChanelLabels = cellstr(num2str((1:size(Winv,1))'));
end

h = figure('units','normalized','outerposition',[0.1 0.1 0.75 0.75],'Visible','off');

CompsToRemove = [];
EEGp = (Winv(:,setdiff(1:size(Winv,2),CompsToRemove))*activations(setdiff(1:size(Winv,2),CompsToRemove),:))';
EEGp = bsxfun(@plus, EEGp, fliplr((1:size(EEGp,2))*nanmean(prctile(EEGp,Offset(2),1)-prctile(EEGp,Offset(1),1))*2 + 0.001));

EEGrange = prctile(EEGp,Offset(2),1)-prctile(EEGp,Offset(1),1);
Yrange = [nanmean(EEGp(:,1)),nanmean(EEGp(:,end))];
Ylim = [Yrange(2)-(nanmean(EEGrange)*2) Yrange(1)+(nanmean(EEGrange)*2)];
ylim(Ylim);

PixSS = get(0,'ScreenSize');
hh = plot(EEGp,'YDataSource','EEGp');
for c = 1:length(hh)
    hh(c).Color = 'k';
end
hhh = gca;
set(hhh,'xtick',1:Srate:size(EEGp,1));
Time = datestr(seconds(0:size(EEGp,1)/Srate),'MM:SS');
set(hhh,'xticklabel',Time)
xlim([1 10*Srate+1]);

ZoomH = zoom(h);
set(ZoomH,'Motion','horizontal');
PanH = pan(h);
set(PanH,'Motion','horizontal');

set(hhh,'ytick',[]);
H = get(hhh,'Children');
highlightChannel(H,flip(ChanelLabels));
linkdata(h);

%% GUI buttons

ButtonPositionsX = 0:round(PixSS(3)*0.75)/4:round(PixSS(3)*0.75);
ButtonPositionsX(1)=ButtonPositionsX(1)+round(0.037*PixSS(4));
ButtonPositionsX(2)=ButtonPositionsX(2)-round(0.037*PixSS(4));
ButtonPositionsX(3)=ButtonPositionsX(3)-round(0.037*PixSS(4));
ButtonPositionsX(4)=ButtonPositionsX(4)-round(0.037*2*PixSS(4));
ButtonPositionsX(5)=ButtonPositionsX(5)-round(0.037*2*PixSS(4));

Buttons(1) = uicontrol('Style', 'pushbutton', 'String', '<<',...
    'Position', [ButtonPositionsX(1) round(0.0185*PixSS(4)) round(0.037*PixSS(4)) round(0.0185*PixSS(4))],...
    'Callback', @move_EEG_backward);

Buttons(2) = uicontrol('Style', 'pushbutton', 'String', '<',...
    'Position', [ButtonPositionsX(2) round(0.0185*PixSS(4)) round(0.037*PixSS(4)) round(0.0185*PixSS(4))],...
    'Callback', @move_EEG_half_backward);

Buttons(4) = uicontrol('Style', 'pushbutton', 'String', '>',...
        'Position', [ButtonPositionsX(4) round(0.0185*PixSS(4)) round(0.037*PixSS(4)) round(0.0185*PixSS(4))],...
        'Callback', @move_EEG_half_forward);

Buttons(5) = uicontrol('Style', 'pushbutton', 'String', '>>',...
        'Position', [ButtonPositionsX(5) round(0.0185*PixSS(4)) round(0.037*PixSS(4)) round(0.0185*PixSS(4))],...
        'Callback', @move_EEG_forward);

Buttons(3) = uicontrol('Style', 'pushbutton', 'String', 'Remove ICs',...
        'Position', [ButtonPositionsX(3)-round(0.222*PixSS(4))/2 round(0.0185*PixSS(4)) round(0.222*PixSS(4)) round(0.0185*PixSS(4))],...
        'Callback', @update_signal);

Buttons(6) = uicontrol('style','edit',...
    'Position', [ButtonPositionsX(3)-round(0.222*PixSS(4))/2 round(0.037*PixSS(4)) round(0.222*PixSS(4)) round(0.0185*PixSS(4))],...
        'Callback', @my_callback);

set(Buttons,'Units','normalized');

function my_callback(h,evt) %#ok<INUSD>
    CompsToRemove = eval(h.String);
end

function update_signal(source, callbackdata) %#ok<INUSD>
    EEGp = (Winv(:,setdiff(1:size(Winv,2),CompsToRemove))*activations(setdiff(1:size(Winv,2),CompsToRemove),:))';
    EEGp = bsxfun(@plus, EEGp, fliplr((1:size(EEGp,2))*nanmean(prctile(EEGp,Offset(2),1)-prctile(EEGp,Offset(1),1))*2 + 0.001));
    EEGrange = prctile(EEGp,Offset(2),1)-prctile(EEGp,Offset(1),1);
    Yrange = [nanmean(EEGp(:,1)),nanmean(EEGp(:,end))];
    Ylim = [Yrange(2)-(nanmean(EEGrange)*2) Yrange(1)+(nanmean(EEGrange)*2)];
    refreshdata(h,'caller');
    msgbox('EEG reconstruction done!')
    ylim(Ylim);
    drawnow;
end

function move_EEG_forward(source, callbackdata) %#ok<INUSD>
    Xlim = xlim;
    CurFrame = round(Xlim(2)-Xlim(1)); % always 30 se
    NextFrame = [Xlim(1)+CurFrame,Xlim(2)+CurFrame];
    xlim(NextFrame);
end

function move_EEG_backward(source, callbackdata) %#ok<INUSD>
    Xlim = xlim;
    CurFrame = round(Xlim(2)-Xlim(1)); % always 30 se
    NextFrame = [Xlim(1)-CurFrame,Xlim(2)-CurFrame];
    xlim(NextFrame);
end

function move_EEG_half_forward(source, callbackdata) %#ok<INUSD>
    Xlim = xlim;
    CurFrame = round((Xlim(2)-Xlim(1))/2); % always 30 se
    NextFrame = [Xlim(1)+CurFrame,Xlim(2)+CurFrame];
    xlim(NextFrame);
end

function move_EEG_half_backward(source, callbackdata) %#ok<INUSD>
    Xlim = xlim;
    CurFrame = round((Xlim(2)-Xlim(1))/2); % always 30 se
    NextFrame = [Xlim(1)-CurFrame,Xlim(2)-CurFrame];
    xlim(NextFrame);
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
%         set(h(a),'userdata',txt{a});
        set(h(a),'linewidth',0.5); % #RM@FBMlab,December 2017
%         set(h(a),'color','k');
        %         bdfcn=['Cp = get(gca,''CurrentPoint''); '  'Xp=Cp(1,1);',  'Yp=Cp(1,2);',   'Zp=Cp(1,3);',  'dat=get(gcbo,''userdata'');',   'ht=text(Xp,Yp,Zp,sprintf(''%s'',dat));',  hTparams]; % #RM@FBMlab,December 2017
        bdfcn=['Cp = get(gca,''CurrentPoint''); '  'Xp=Cp(1,1);',  'Yp=Cp(1,2);',   'Zp=Cp(1,3);',  'dat=get(gcbo,''userdata'');',   'ht=text(Xp,Yp,Zp,sprintf(''%s'',dat)); if get(gco,''linewidth'')==0.5,set(gco,''linewidth'',3),set(gco,''color'',''red''),else set(gco,''linewidth'',0.5),set(gco,''color'',''k''),end; ',  hTparams]; % #RM@FBMlab,December 2017
        set(h(a),'buttondownfcn',bdfcn);
    end
else
%     set(h,'userdata',txt);
    set(h,'linewidth',0.5); % #RM@FBMlab,December 2017
%     set(h,'color','k');
    %     bdfcn=['Cp = get(gca,''CurrentPoint''); '  'Xp=Cp(1,1);',  'Yp=Cp(1,2);',   'Zp=Cp(1,3);', 'dat=get(gcbo,''userdata'');', 'ht=text(Xp,Yp,Zp,sprintf(''%s'',dat));',  hTparams]; % #RM@FBMlab,December 2017
    bdfcn=['Cp = get(gca,''CurrentPoint''); '  'Xp=Cp(1,1);',  'Yp=Cp(1,2);',   'Zp=Cp(1,3);', 'dat=get(gcbo,''userdata'');', 'ht=text(Xp,Yp,Zp,sprintf(''%s'',dat)); if get(gco,''linewidth'')==0.5,set(gco,''linewidth'',3),set(gco,''color'',''red''),else set(gco,''linewidth'',0.5),set(gco,''color'',''k''),end; ',  hTparams]; % #RM@FBMlab,December 2017
    set(h,'buttondownfcn',bdfcn);
end

end
