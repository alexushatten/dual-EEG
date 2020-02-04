function [Shifts,ToCheck] = spike_aligner_v3( EEG, Marker_Time, ChannelToShow, ChannelLabel, FrameLength )
% spike_aligner_v2: align epochs of spikes with reference to channel with
%   maximal absolute amplitude
%
% [Shifts,ToCheck] = spike_aligner_v3( EEG, Marker_Time,...
%                    ChannelLabel, ChannelToShow, FrameLength )
%
%  Inputs
% --------
% EEG : [channels x time] EEG array
% Marker_Time : [n x 1] subscripts indices (time frames) of markers to adjust
% ChannelToShow : [1 x 1] subscript index of channel to show
% ChannelLabel : char, label of channel to show
% FrameLength : numeric, length of frame to display (in samples)
%
%  Outputs
% ---------
% Shifts : shifts of signal (in samples) for each epoch
% ToCheck: logical, epochs marked during alignment as "to be checked later"
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, April 2019, updated May 2019
%-------------------------------------------------------------------------

% loop across epochs
Shifts = nan(length(Marker_Time),1);
ToCheck = false(length(Marker_Time),1);
Exit = false;
m = 1;
while m<=length(Marker_Time) && ~Exit
    Frames = (Marker_Time(m)-round(FrameLength/2)):(Marker_Time(m)+round(FrameLength/2));
    if Exit
        break;
    end
%     RefX = repmat(round(FrameLength/2)+1,FrameLength+1,1);
    RefY = linspace(min(vect(EEG(ChannelToShow,Frames))),max(vect(EEG(ChannelToShow,Frames))),FrameLength+1);
    
    [Shifts(m), Exit, ToCheck(m), Prev] = overlay_spike_on_avg_spike(EEG(ChannelToShow,Frames)',RefY*5,ChannelLabel,ToCheck(m),m,length(Marker_Time));
    if Prev
        if m == 1
            asdf = msgbox('First epoch, cannot go further back !');
            waitfor(asdf);
        else
            m = m-1;
        end
    else
        m = m+1;
    end
end

end

function [Shift, Exit, ToCheck, Prev] = overlay_spike_on_avg_spike(EEG,AvgSpike,ChanelLabel,ToCheck,Current,Total)

Prev = false;
Exit = false;
% make figure
h = figure('units','normalized','outerposition',[0.1 0.1 0.75 0.75],'name','Spike aligner v3',...
    'Visible','on');
hold on;
% PixSS = get(0,'ScreenSize'); %PixSS(3) = PixSS(3)/1.8750;

%% plot
RefAvg = [nan(1,size(AvgSpike,2)),AvgSpike,nan(1,size(AvgSpike,2))];
hh = plot(repmat(size(EEG,1)+round(size(EEG,1)/2),length(RefAvg),1),RefAvg,'color',[1 1 1]);
title([ChanelLabel,' (polarity upward)', [num2str(Current), ' / ', num2str(Total)]]);
xlabel({'Use arrow keys to shift / scale current spike,',...
    '[Enter] to continue',...
    '[Escape] to exit spike alignment',...
    '[C] for marking as "to check"',...
    '[P] for polarity switch',...
    '[B] for going back to previous spike',...
    'and [R] for reset current alignment'})
G = ylabel('OK');
set(G,'rotation',0);

RefLineXavg = [size(EEG,1)+round(size(EEG,1)/2) size(EEG,1)+round(size(EEG,1)/2)];
RefLineXcurrent = RefLineXavg;
Ylim = ylim;
MaxVal = max(AvgSpike)*20;
line(RefLineXavg,[-MaxVal MaxVal],'color','red');
hRefLineXcurrent = line(RefLineXcurrent,[-MaxVal MaxVal],'color','blue');

Pswitch = false;
hh.YDataSource = 'RefAvg';
Y = [nan(size(EEG,1),1);EEG;nan(size(EEG,1),1)];
Y_BKP = Y;
hhh = plot(Y,'k');
hhh.YDataSource = 'Y';
Y = [nan(size(EEG,1),1);EEG;nan(size(EEG,1),1)];

refreshdata(hhh,'caller');
xlim([size(EEG,1) size(EEG,1)*2]);
ylim(Ylim);
set(gca,'xtick',0:100:size(EEG,1)*3);
set(gca,'xticklabel',{-size(EEG,1)+1:100:size(EEG,1)*2});
drawnow;

%% init
Shift = 0;

%% GUI button
% h.Visible = 'off';
% ButtonPositionsX = 0:round(PixSS(3)*0.75)/3:round(PixSS(3)*0.75);
% ButtonPositionsX(1)=ButtonPositionsX(1)+round(0.037*PixSS(4));
% ButtonPositionsX(2)=ButtonPositionsX(2)-round(0.037*PixSS(4));
% ButtonPositionsX(3)=ButtonPositionsX(3)-round(0.037*PixSS(4));
% ButtonPositionsX(4)=ButtonPositionsX(4)-round(0.037*4*PixSS(4));

% Buttons(1) = uicontrol('Style', 'pushbutton', 'String', 'Check later!',...
%     'Position', [ButtonPositionsX(1) round(0.0185*PixSS(4)) round(1/9*PixSS(4)) round(0.0185*PixSS(4))],...
%     'Callback', @check_it_later);
% 
% Buttons(2) = uicontrol('Style', 'pushbutton', 'String', 'Reset current',...
%     'Position', [ButtonPositionsX(4) round(0.0185*PixSS(4)) round(1/9*PixSS(4)) round(0.0185*PixSS(4))],...
%     'Callback', @reset_shift);

% set(Buttons,'Units','normalized');

%% keyboard controls
set(gcf,'KeyPressFcn',@ReadKey);

function check_it_later(source, callbackdata) %#ok<INUSD>
    if ToCheck == true;
        ToCheck = false;
        set(G,'color',[0.15 0.15 0.15]);
        set(G,'string','OK');
        set(G,'fontweight','normal');
    else
        ToCheck = true;
        set(G,'color',[1 0 0]);
        set(G,'string','TO CHECK');
        set(G,'fontweight','bold');
    end
end

function reset_shift(source, callback) %#ok<INUSD>
    RefLineXcurrent = RefLineXavg;
    Y = Y_BKP;
    delete(hRefLineXcurrent); refresh(h);
    hRefLineXcurrent = line(RefLineXcurrent,[-MaxVal MaxVal],'color','blue');
    refreshdata(hhh,'caller');
    Shift = 0;
end

function ReadKey(source, callback) %#ok<INUSD>

Key=get(gcf,'CurrentCharacter');

switch double(Key)
    case 28 % left arrow => move left ...
        % => pan at the end and remove values at the beginning
        Y(end+1) = nan;
        Y(1) = [];
        RefLineXcurrent = RefLineXcurrent-1;
        delete(hRefLineXcurrent); refresh(h);
        hRefLineXcurrent = line(RefLineXcurrent,[-MaxVal MaxVal],'color','blue');
        Shift = Shift-1;
        refreshdata(hhh,'caller');
    case 29 % right arrow => move right ...
        % => pan at the beginning and remove values at the end
        Temp = [nan;Y];
        Temp(end) = [];
        Y = Temp;
        RefLineXcurrent = RefLineXcurrent+1;
        delete(hRefLineXcurrent); refresh(h);
        hRefLineXcurrent = line(RefLineXcurrent,[-MaxVal MaxVal],'color','blue');
        Shift = Shift+1;
        refreshdata(hhh,'caller');
    case 13 % enter, continue with next spike
        close(h);
    case 27 % escape, exit spike alignment
        close(h);
        Exit = true;
        return;
    case 30 % arrow up, scale up current spike
        Y = Y*1.05;
        refreshdata(hhh,'caller');
    case 31 % arrow down, scale down current spike
        Y = Y*0.95;
        refreshdata(hhh,'caller');
    case 112 % key "P", switch polarity
        if Pswitch
            Pswitch = false;
        else
            Pswitch = true;
        end
        if Pswitch
            ParseChanLab = regexp(ChanelLabel,'-','split');
            if length(ParseChanLab)>1
                NewChanLab = cellstr([ParseChanLab{1}{2},'-',ParseChanLab{1}{1}]);
                title([NewChanLab,' (polarity switch on current spike)',[num2str(Current), ' / ', num2str(Total)]]);
            else
                title([ChanelLabel,' (polarity switch on current spike)',[num2str(Current), ' / ', num2str(Total)]]);
            end
        else
            title([ChanelLabel,' (default polarity)',[num2str(Current), ' / ', num2str(Total)]]);
        end
        Y = -Y;
        refreshdata(hhh,'caller');
    case 114 % key "R", reset alignment
        reset_shift;
    case 99 % key "C", check spike later!
        check_it_later;
    case 98 % key "B", go back to previous spike epoch
        Prev = true;
        close(h);
        return;
end
end

%% show figure
% h.Visible = 'on';
waitfor(h);

end
