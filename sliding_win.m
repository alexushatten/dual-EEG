function Out = sliding_win( Data, Window , Fun, Graph )
%sliding_win: backward and forward sliding window
%
% Usage:
%-------------------------------------------------------------------------
% 
% Out = sliding_win( Data, Window , Fun, Graph )
%
% Inputs:
%-------------------------------------------------------------------------
%
% Data : signal, n x 1 array
%
% Window : window size in samples
%
% Fun: function handle, e.g. @mean
%
% Graph: plot graph or not
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, March 2018
%-------------------------------------------------------------------------

if isempty(Graph)==1
    Graph = 1;
end

LengthData = length(Data);

if all(~isinf(Data))==0 % Check for Infinite values in the data (typically present in some reaction times data when the subject didn't respond)
    Data(isinf(Data))=nan(1,1); % Transform it in NaN values
end

for i = 1+round(Window/2):LengthData-(round(Window/2)-1)
    Period(:,i-(Window/2)) = Data(i-round(Window/2):i+(round(Window/2)-1));
end

% Smoothed = mean(Period);
% Download the function ignoreNaN on Matlab File Exchange and use this
% instead:

Out = ignoreNaN(Period,Fun)'; % Ignore NaN values to perform the mean

if Graph == 1
    figure;
    plot(1:size(Data,2),Data);
    hold on;
    plot(round(Window/2):(size(Out,2)+round(Window/2)-1),Out,'r');
end

end

