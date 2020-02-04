function remove_duplicates_mrk_file( MrkFile )
% remove_duplicates_mrk_file: removes duplicate markers from specified
% file (i.e. same onset, offset and label)
%
% remove_duplicates_mrk_file( MrkFile )
%
%  Inputs
% --------
% MrkFile: .mrk file path
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, October 2018
%-------------------------------------------------------------------------

try
    [tS, tE, tL] = read_mrk_file(MrkFile);
catch ME
    warning([ME.message,' Trying fallback marker reading function instead.'])
    [tS, tE, tL] = read_mrk_Cartool(MrkFile);
end

JokerSet = {'_','$','!','?','#','/'};
% Find Joker
js = 1;
Limit = length(JokerSet);
while any(~cellfun(@isempty,regexp(tL,JokerSet{js}))) && js<Limit
    js = js + 1;
end
if js==Limit
    % start combining
    NumChar = 2;
    
    JokerSetChar = char(JokerSet)';
    CombJoker = cellstr(JokerSetChar(perms(1:NumChar)));    
    js = 1;
    while any(~cellfun(@isempty,regexp(tL,CombJoker{js}))) && ~(NumChar>length(JokerSet))
        Limit = length(CombJoker);
        while any(~cellfun(@isempty,regexp(tL,CombJoker{js}))) && js<Limit
            js = js + 1;
        end
        NumChar = NumChar + 1;
        if NumChar>length(JokerSet)
            break;
        end
        CombJoker = cellstr(JokerSetChar(perms(1:NumChar))); 
    end
    if js==Limit && NumChar>length(JokerSet)
        error('Could not find set of characters not already present in the marker labels in order to decombine them!')
    else
        Joker = CombJoker{js};
    end
else
    Joker = JokerSet{js};
end 

CombMrk = cellstr([num2str(tS),repmat(Joker,length(tL),1),num2str(tE),repmat(Joker,length(tL),1),char(tL)]);

CheckCombMrk = match_vectors(CombMrk,CombMrk,1);
if isa(CheckCombMrk,'cell')
    CombMrk = unique(CombMrk);
    
    DecombMrk = regexp(CombMrk,Joker,'split');
    [tSnew,tEnew] = deal(nan(size(DecombMrk)));
    for m = 1:length(DecombMrk)
        tSnew(m) = str2num(DecombMrk{m}{1}); %#ok<ST2NM>
        tEnew(m) = str2num(DecombMrk{m}{2}); %#ok<ST2NM>
        tLnew{m} = DecombMrk{m}{3}; %#ok<AGROW>
    end
    tLnew = tLnew';
    copyfile(MrkFile,spm_file(MrkFile,'suffix','_with_duplicates'));
    write_mrk_file_Cartool(MrkFile,tSnew,tEnew,tLnew);
else
    fprintf('No duplicate markers found!\n')
end

end

