function [ Hits, FP, Misses, YinT, TinY ] = timings_matching( T, Y, EOF, Lag )
% timings_matching: counts number of events that co-occur in two vectors,
% with a certain sample lag allowed
%
% [ Hits, FP, Misses, YinT, TinY] = timings_matching( T, Y, EOF, Lag )
%
%  Inputs
% --------
% T  : nx1 vector, ground truth, vector of timings, if nx2 the first column
%       is assumed to represent the start and the second the end of the event.
%
% Y  : nx1 vector, vector of timings to test, if nx2 the first column is
%       assumed to represent the start and the second the end of the event.
%
% EOF : integer, end of signals (normally both T and Y should end at same time!).
%
% Lag : allowed lag (in samples) between timings of T and Y.
%
%
%  Outputs
% ---------
% Hits   : number of true positives
%
% FP     : number of false positives
%
% Misses : number of false negatives
%
% YinT   : timings of Y that match in T (if NaN it is not found in T)
%
% TinY   : timings of T that match in Y (if NaN it is not found in Y)
%
%
%-------------------------------------------------------------------------
% Nota bene:
%           Hits + FP = length(Y)
%
%  However, depending on Lag value, some events might be merged...
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, June 2018
%-------------------------------------------------------------------------

% might happen...
T = T(~(sum(isnan(T),2)>0),:);
Y = Y(~(sum(isnan(Y),2)>0),:);

if isempty(T)
    Hits = nan;
    FP = nan;
    Misses = nan;
    YinT = nan;
    TinY = nan;
    warning('No timings in reference vector, exiting...');
    return;
elseif isempty(Y)
    Hits = nan;
    FP = nan;
    Misses = nan;
    YinT = nan;
    TinY = nan;
    warning('No timings in test vector, exiting...');
    return;
end

if size(T,2)==2 && size(Y,2)==2
    DurationFlag = true;
elseif size(T,2)==1 && size(Y,2)==1
    DurationFlag = false;
else
    error('Y1 and Y2 do not have proper dimensions, both should be n x 1 or n x 2')
end

% EOF = max([Y1eof;Y2eof]);
if DurationFlag
    Y1sTemp = false(EOF,1);
    for t = 1:size(T,1)
        Y1sTemp(T(t,1):T(t,2))=true;
    end
    Y1sTemp = find(Y1sTemp);
    Y2sTemp = false(EOF,1);
    for t = 1:size(Y,1)
        Y2sTemp(Y(t,1):Y(t,2))=true;
    end
    Y2sTemp = find(Y2sTemp);
else
    Y1sTemp = T;
    Y2sTemp = Y;
end

Y1s = smooth(single(dnif(Y1sTemp,EOF)),Lag*2)>0;
Y2s = smooth(single(dnif(Y2sTemp,EOF)),Lag*2)>0;

% S = bin2dec([num2str(Y1s),num2str(Y2s)]);

Onsets1 = find(diff(Y1s)>0);
Offsets1 = find(diff(Y1s)<0);
Onsets2 = find(diff(Y2s)>0);
Offsets2 = find(diff(Y2s)<0);

if length(Onsets1)>length(Offsets1)
    Offsets1(end+1) = EOF;
end
if length(Onsets2)>length(Offsets2)
    Offsets2(end+1) = EOF;
end
if length(Offsets1)>length(Onsets1)
    Onsets1BKP = Onsets1;
    Onsets1 = nan(size(Onsets1)+[1 0]);
    Onsets1(2:end) = Onsets1BKP;
    Onsets1(1) = 0;
end
if length(Offsets2)>length(Onsets2)
    Onsets2BKP = Onsets2;
    Onsets2 = nan(size(Onsets2)+[1 0]);
    Onsets2(2:end) = Onsets2BKP;
    Onsets2(1) = 0;
end

Found1 = false(length(Onsets1),1);
TinY = nan(length(Onsets1),1);
for t = 1:length(Onsets1)
    Temp1 = Onsets1(t)>=Onsets2;
    Temp2 = Offsets1(t)<=Offsets2;
    Temp3 = Onsets1(t)<=Offsets2;
    Temp4 = Offsets1(t)>=Onsets2;
    if any((Temp1+Temp2)>1) || any((Temp3+Temp4)>1);
        Where1temp = unique([find((Temp1+Temp2)>1);find((Temp3+Temp4)>1)]);
        TinY(t) = Onsets2(Where1temp(~isempty(Where1temp)))+Offsets2(Where1temp(~isempty(Where1temp)));
        Found1(t) = true;
    end
end

Found2 = false(length(Onsets2),1);
YinT = nan(length(Onsets2),1);
for t = 1:length(Onsets2)
    Temp1 = Onsets2(t)>=Onsets1;
    Temp2 = Offsets2(t)<=Offsets1;
    Temp3 = Onsets2(t)<=Offsets1;
    Temp4 = Offsets2(t)>=Onsets1;
    if any((Temp1+Temp2)>1) || any((Temp3+Temp4)>1);
        Where2temp = unique([find((Temp1+Temp2)>1);find((Temp3+Temp4)>1)]);
        YinT(t) = Onsets1(Where2temp(~isempty(Where2temp)))+Offsets1(Where2temp(~isempty(Where2temp)));
        Found2(t) = true;
    end
end

% Y1hits = sum(Found1);
Misses = sum(~Found1);
Hits = sum(Found2);
FP = sum(~Found2);

end

