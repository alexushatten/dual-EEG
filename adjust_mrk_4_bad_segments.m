function [NewMrkOnsets, NewMrkOffsets, NewMrkLabels] = adjust_mrk_4_bad_segments( MrkOnsets, MrkOffsets, MrkLabels, BadOnsets, BadOffsets )
% adjust_mrk_4_bad_segments: adjust marker timings when excising bad
% segments from recording, i.e. shift markers after bad segments, exclude
% markers inside bad segments, and adjust onsets, offsets and duration of
% markers partially overlapping with bad segments
%
% NewMarkers = adjust_mrk_4_bad_segments( Onsets, Offsets, Labels, BadOnsets, BadOffsets )
%
% NewMarkers = adjust_mrk_4_bad_segments( Onsets, Offsets, Labels, Bad_Segments )
%
%  Inputs
% --------
% MrkOnsets: [m x 1] integer, onset of markers
% MrkOffsets: [m x 1] integer, offset of markers
% MrkLabels: {m x 1} cell of strings, label of markers
%
% BadOnsets: [b x 1] integer, onset of bad segments
% BadOffsets: [b x 1] integer, offset of bad segments
% 
% Alternatively to specifying onsets and offsets of bad segments, you can
% provide a logical vector:
% Bad_Segments: [L x y] logical, whether each time frame is bad (1) or good (0)
%
%  Outputs
% ---------
% NewMrkOnsets: [n x 1] integer, adjusted onset of markers
% NewMrkOffsets: [n x 1] integer, adjusted offset of markers
% NewMrkLabels: {n x 1} cell of strings, adjusted labels of markers
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, May 2019
%-------------------------------------------------------------------------

if islogical(BadOnsets)
    BadIsLogical = true;
    Bad = BadOnsets;
    % Construct onsets and offsets of bad segments:
    BadOnsets = find(diff(Bad)>0);
    BadOffsets = find(diff(Bad)<0);
    % Correct for BOF & EOF:
    if length(BadOnsets)>length(BadOffsets)
        BadOffsets(end+1) = size(Bad,2);
    end
    if length(BadOffsets)>length(BadOnsets)
        BadOnsetsBKP = BadOnsets;
        BadOnsets = nan(size(BadOnsets)+[1 0]);
        BadOnsets(2:end) = BadOnsetsBKP;
        BadOnsets(1) = 1;
    end
    if Bad(1) && Bad(end) % in case both BOF & EOF are bad...
        BadOffsets(end+1) = size(Bad,2);
        
        BadOnsetsBKP = BadOnsets;
        BadOnsets = nan(size(BadOnsets)+[1 0]);
        BadOnsets(2:end) = BadOnsetsBKP;
        BadOnsets(1) = 1;
    end
end

if length(MrkOnsets)~=length(MrkOffsets) || length(MrkOnsets)~=length(MrkLabels) || length(BadOnsets)~=length(BadOffsets)
    error('Input size mismatch')
end

BadTemp = [];
Count = 0;
for b = 1:length(BadOnsets)
    Ltemp = BadOffsets(b)-BadOnsets(b)+1;
    BadTemp(Count+1:Count+Ltemp) = BadOnsets(b):BadOffsets(b);
    Count = Count+Ltemp;
end
L = max(MrkOffsets(end),BadOffsets(end));
Bad = dnif(BadTemp,L);

NewMrkOnsets = [];
NewMrkOffsets = [];
NewMrkLabels = {};
Count = 1;
for m = 1:length(MrkOnsets)
    % Develop marker:
    MrkTemp = MrkOnsets(m):MrkOffsets(m);
    Remains = setdiff(MrkTemp,BadTemp);
    if ~isempty(Remains)
        BadBefore = sum(Bad(1:Remains));
        NewMrkTemp = Remains-BadBefore;
        NewMrkOnsets(Count) = NewMrkTemp(1); %#ok<AGROW>
        NewMrkOffsets(Count) = NewMrkTemp(end); %#ok<AGROW>
        NewMrkLabels(Count) = MrkLabels(m); %#ok<AGROW>
        Count = Count+1;
    end
end

NewMrkOnsets = NewMrkOnsets(:);
NewMrkOffsets = NewMrkOffsets(:);
NewMrkLabels = NewMrkLabels(:);

end

% BadOnsets_BKP = BadOnsets;
% BadOffsets_BKP = BadOffsets;
% NewMrkOnsets = MrkOnsets;
% NewMrkOffsets = MrkOffsets;
% NewMrkLabels = MrkLabels;
%
% % The following code looks cumbersome but is actually the only way to
% % clearly disentangle each situation, because the two only things we can
% % assume are the following:
% %   - The onset of a marker is always after the offset of the same marker;
% %   - Bad segments never overlap with themselve.
% % Then there are corollaries, but looking at overlap using e.g. function
% % "dnif" will not work because for markers with duration > 0 we
% % need to reconstruct from scratch and there can be e.g. multiple markers
% % at the same time, with the same duration...
%
% for b = 1:length(BadOnsets)
%
%     Cases = sum([(NewMrkOnsets>=BadOnsets(b)),...
%         2*(NewMrkOnsets<=BadOffsets(b)),...
%         4*(NewMrkOffsets>=BadOnsets(b)),...
%         8*(NewMrkOffsets<=BadOffsets(b))],2);
%
%     % Responses: SUM: meaning:
%     % [0 0 0 0] => 0: mrk onset before bad onset and after bad offset (impossible)
%     % [1 0 0 0] => 1: mrk onset after bad onset and after bad offset (after?), mrk offset before bad onset and after bad offset => (impossible)
%     % [0 2 0 0] => 2: mrk onset before bad onset and before bad offset (before?), mrk offset before bad onset and after bad offset => (impossible)
%     % [1 2 0 0] => 3: mrk onset after bad onset and before bad offset (middle?), mrk offset before bad onset and after bad offset => (impossible)
%     %
%     % [0 0 4 0] => 4: mrk onset before bad onset and after bad offset (impossible)
%     % [1 0 4 0] => 5: mrk onset after bad onset and after bad offset (after?), mrk offset after bad onset and after bad offset => CASE "AFTER"!
%     % [0 2 4 0] => 6: mrk onset before bad onset and before bad offset (before?), mrk offset after bad onset and after bad offset => CASE "SPAN"!
%     % [1 2 4 0] => 7: mrk onset after bad onset and before bad offset (middle?), mrk offset after bad onset and after bad offset => CASE "END-OVERLAP"!
%     %
%     % [0 0 0 8] => 8: mrk onset before bad onset and after bad offset (impossible)
%     % [1 0 0 8] => 9: mrk onset after bad onset and after bad offset (after?), mrk offset before bad onset and before bad offset => (impossible)
%     % [0 2 0 8] => 10: mrk onset before bad onset and before bad offset (before?), mrk offset before bad onset and before bad offset => CASE "BEFORE"!
%     % [1 2 0 8] => 11: mrk onset after bad onset and before bad offset (middle?), mrk offset before bad onset and before bad offset => (impossible)
%     %
%     % [0 0 4 8] => 12: mrk onset before bad onset and after bad offset (impossible)
%     % [1 0 4 8] => 13: mrk onset after bad onset and after bad offset (after?), mrk offset after bad onset and before bad offset => (impossible)
%     % [0 2 4 8] => 14: mrk onset before bad onset and before bad offset (before?), mrk offset after bad onset and before bad offset => CASE "START-OVERLAP"!
%     % [1 2 4 8] => 15: mrk onset after bad onset and before bad offset (middle?), mrk offset after bad onset and before bad offset => CASE "INSIDE"!
%
%     % 5 => "AFTER" => subtract duration of bad +1
%     %                 from start and end
%     % 6 => "SPAN" => start ok, subtract duration
%     %                 of bad +1 from end
%     % 7 => "END-OVERLAP" => start ok, end becomes
%     %                 start of bad -1
%     % 10 => "BEFORE" => start and end ok!
%     % 14 => "START-OVERLAP" => start becomes start
%     %                 of bad -1, subtract duration
%     %                 of bad +1 from end
%     % 15 => "INSIDE" => NaN them, remove at the end!
%
%     if ~isempty(intersect([1,2,3,4,8,9,11,12,13],Cases))
%         error('Check input markers & bad segments onsets / offsets, some offsets seem to lie before the onsets...')
%     elseif ~all(isnan(NewMrkOnsets(Cases==0))) || ~all(isnan(NewMrkOffsets(Cases==0)))
%         % Case 0 can occur when markers inside bad segments are NaNed, but
%         % otherwise it should never occur !
%         error('Check input markers & bad segments onsets / offsets, some offsets seem to lie before the onsets...')
%     end
%
%     NewMrkOnsets(Cases==5) = NewMrkOnsets(Cases==5)-(BadOffsets(b)-BadOnsets(b)+1);
%     NewMrkOffsets(Cases==5) = NewMrkOffsets(Cases==5)-(BadOffsets(b)-BadOnsets(b)+1);
%
% %     NewMrkOnsets(Cases==6) % ok, nothing to do
%     NewMrkOffsets(Cases==6) = NewMrkOffsets(Cases==6)-(BadOffsets(b)-BadOnsets(b)+1);
%
% %     NewMrkOnsets(Cases==7) % ok, nothing to do
%     NewMrkOffsets(Cases==7) = BadOnsets(b)-1;
%
% %     NewMrkOnsets(Cases==10) % ok, nothing to do
% %     NewMrkOffsets(Cases==10) % ok, nothing to do
%
%     NewMrkOnsets(Cases==14) = BadOnsets(b)-1;
%     NewMrkOffsets(Cases==14) = NewMrkOffsets(Cases==14)-(BadOffsets(b)-BadOnsets(b)+1);
%
%     NewMrkOnsets(Cases==15) = nan;
%     NewMrkOffsets(Cases==15) = nan;
%
%     % Update bad segments' onsets and offsets
%     BadOnsets = BadOnsets-(BadOffsets(b)-BadOnsets(b)+1);
%     BadOffsets = BadOffsets-(BadOffsets(b)-BadOnsets(b)+1);
% end
%
% ToRemove = isnan(NewMrkOnsets);
% NewMrkOnsets(ToRemove) = [];
% NewMrkOffsets(ToRemove) = [];
% NewMrkLabels(ToRemove) = [];
