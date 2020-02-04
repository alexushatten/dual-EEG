function [Shifts,ToCheck] = spike_aligner_hdEEG( File, MarkerOfInterest, OldFile, OldMarkerOfInterest, ChannelOfInterest, FrameLength, Banana )
% spike_aligner_hdEEG: align epochs of spikes with reference to channel with
%   maximal absolute amplitude
%
% [Shifts, ToCheck] = spike_aligner_hdEEG( File, MarkerOfInterest, ChannelOfInterest, FrameLength )
%
%  Inputs
% --------
% File :                char, path to .sef file
% MarkerOfInterest :    char, markers to consider for alignment (other markers
%                       are ignored)
% OldFile:              char, path to reference .sef file for which marking
%                       was done (useful for quick and dirty averaging)
% OldMarkerOfInterest:  char, how the marker of interest was labeled before
%                       new realignment strategy
% ChannelOfInterest :   char, channel name to consider for alignment (other
%                       channels are not displayed)
% FrameLength :         numeric, length of frame to display (in samples)
% Banana: string, either 'double' or 'triple'
%
%  Outputs
% ---------
% Shifts : shifts of signal (in samples) for each epoch
% ToCheck: logical, epochs marked during alignment as "to be checked later"
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, April 2019
%-------------------------------------------------------------------------

% Get EEG traces
[EEG, Hdr] = dual_load_sef(File);
[~,labelsB,~,~,EEGb] = bipolar_montage_scalp(Banana,EEG);
% quick filter + demean (for intracranial it is not necessary but for scalp it is!)
EEGb = dual_bw_bp_filt( EEGb, Hdr.samplingfreq, 1, 70 );
EEGb = detrend(brainstorm_notch(EEGb,Hdr.samplingfreq,50)')';

%% for old average plot

if ~strcmp(File,OldFile) % in case we just want the average out of the same EEG file and not another
    [asdf,asdf_hdr] = dual_load_sef(OldFile);
    [~,bip_labels_m,~,~,asdfb] = bipolar_montage_scalp(Banana,asdf);
else
    asdfb = EEGb; asdf_hdr = Hdr;
end
[oldMarker_T1,~,oldMarker_Label] = read_mrk_Cartool([OldFile,'.mrk']);
oldMarkersToDo = strcmp({OldMarkerOfInterest},oldMarker_Label);
AllMatches = oldMarker_T1(oldMarkersToDo);
if isempty(AllMatches)
    warning('No matching marker was found, please pick another file or change marker label');
    Shifts = nan; ToCheck = nan;
    return;
end

Epochs = nan(size(asdfb,1),round(FrameLength),length(AllMatches));
% quick filter + demean
asdfb = dual_bw_bp_filt( asdfb, asdf_hdr.samplingfreq, 1, 70 );
asdfb = detrend(brainstorm_notch(asdfb,asdf_hdr.samplingfreq,50)')';
for ep = 1:length(AllMatches)
    try
        Epochs(:,:,ep) = asdfb(:,AllMatches(ep)-round(FrameLength/2):AllMatches(ep)+round(FrameLength/2)-1);
    catch
        % in case there are some spikes too close to BOF / EOF (beginning
        % or end of recording... but cannot be both anyhow)
        if (AllMatches(ep)-round(FrameLength/2))<1 % too close to BOF
            TempEpoch = asdfb(:,1:AllMatches(ep)+round(FrameLength/2)-1);
            Epochs(:,:,ep) = [nan(size(asdfb,1),abs(AllMatches(ep)-round(FrameLength/2))+1),TempEpoch];
%             Epochs(:,:,ep) = [nan(size(asdfb,1),(AllMatches(ep)-round(FrameLength/2))-1),TempEpoch];
        end
        if (AllMatches(ep)+round(FrameLength/2)-1)>size(asdfb,2) % too close to EOF
            TempEpoch = asdfb(:,AllMatches(ep)-round(FrameLength/2):end);
            Epochs(:,:,ep) = [TempEpoch,nan(size(asdfb,1),(AllMatches(ep)+round(FrameLength/2)-1)-size(asdfb,2))];
        end
    end    
end
% quick and dirty averaging
% AvgSpike = mean(Epochs,3);
AvgSpike = nanmean(Epochs,3); % in case there are some spikes too close to BOF / EOF (beginning or end of recording)

% [~,asdfidx] = max(max(AvgSpike,[],2));

%% quickly checking whether the channel selected is appropriate

ChannelToShow = strcmp({ChannelOfInterest},bip_labels_m');

h = plot_spike( AvgSpike', find(ChannelToShow), asdf_hdr.samplingfreq, bip_labels_m' );
title(regexprep(OldMarkerOfInterest,'_',' '));
waitfor(h);

% % The following cannot be executed because dat is accessible in the main
% workspace, not the workspace of the function
% if ~isempty(dat)
%     if ~strcmp(ChannelOfInterest,dat)
%         warning('Changing channel to the one selected')
%         ChannelOfInterest = dat;
%     end
% end

%% for new timings

% Matching channel
ChannelToShow = strcmp({ChannelOfInterest},labelsB');

% Get markers
[Marker_T1,Marker_T2,Marker_Label] = read_mrk_Cartool([File,'.mrk']);

% Matching markers and get center between start and end
MarkersToDo = strcmp({MarkerOfInterest},Marker_Label);
if sum(MarkersToDo)==0
    warning('No marker label matching the one of interest was found, exiting...');
    Shifts = nan;
    ToCheck = nan;
    return;
end
Marker_Time = round((Marker_T1(MarkersToDo)+Marker_T2(MarkersToDo))/2);

% loop across epochs
Shifts = nan(length(Marker_Time),1);
ToCheck = false(length(Marker_Time),1);
Exit = false;
m = 1;
ScaleFactor = 1;
while m<=length(Marker_Time) && ~Exit
    Frames = (Marker_Time(m)-round(FrameLength/2)):(Marker_Time(m)+round(FrameLength/2));
    if Exit
        break;
    end
%     RefX = repmat(round(FrameLength/2)+1,FrameLength+1,1);
%     RefY = linspace(min(vect(EEG(ChannelToShow,Frames))),max(vect(EEG(ChannelToShow,Frames))),FrameLength+1);
    
    try
        DisplayMe = EEGb(ChannelToShow,Frames);
    catch
        % in case there are some spikes too close to BOF / EOF (beginning
        % or end of recording... but cannot be both anyhow)
        if any(Frames<1) % too close to BOF
            TempEpoch = EEGb(ChannelToShow,1:Frames(end));
            DisplayMe = [nan(1,abs(Marker_Time(m)-round(FrameLength/2))),TempEpoch];
        end
        if any(Frames>size(EEGb,2)) % too close to EOF
            TempEpoch = EEGb(ChannelToShow,Frames(1):end);
            DisplayMe = [TempEpoch,nan(1,(Marker_Time(m)+round(FrameLength/2))-size(EEGb,2)-1)];
        end
    end

    [Shifts(m), Exit, ToCheck(m), Prev, ScaleFactor] = overlay_spike_on_avg_spike(DisplayMe',AvgSpike(ChannelToShow,:),labelsB(ChannelToShow),ToCheck(m),m,length(Marker_Time),File,MarkerOfInterest, ScaleFactor);
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


function [Shift, Exit, ToCheck, Prev, ScaleFactor] = overlay_spike_on_avg_spike(EEG,AvgSpike,ChanelLabel,ToCheck,Current,Total,File,MarkerOfInterest, ScaleFactor)

Prev = false;
Exit = false;
% make figure
% h = figure('units','normalized','outerposition',[0.1 0.1 0.75 0.75],...
%     'Visible','on');
% PixSS = get(0,'ScreenSize'); %PixSS(3) = PixSS(3)/1.8750;

h = figure('units','normalized','outerposition',[0.1 0.1 0.75 0.75],'name',File,...
    'Visible','on');
hold on;

%% plot
RefAvg = [nan(1,size(AvgSpike,2)),AvgSpike,nan(1,size(AvgSpike,2))];
hh = plot(RefAvg,'color',[0.5 0.5 0.5]);
title([ChanelLabel,' (polarity upward)', regexprep(MarkerOfInterest,'_|^',' '), [num2str(Current), ' / ', num2str(Total)]]);
% title([ChanelLabel,' (polarity upward)', [num2str(Current), ' / ', num2str(Total)]]);
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
Y = [nan(size(EEG,1),1);EEG*ScaleFactor;nan(size(EEG,1),1)];
Y_BKP = Y;
hhh = plot(Y,'k');
hhh.YDataSource = 'Y';
Y = [nan(size(EEG,1),1);EEG*ScaleFactor;nan(size(EEG,1),1)];

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
        ScaleFactor = ScaleFactor*1.05;
        refreshdata(hhh,'caller');
    case 31 % arrow down, scale down current spike
        Y = Y*0.95;
        ScaleFactor = ScaleFactor*0.95;
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
