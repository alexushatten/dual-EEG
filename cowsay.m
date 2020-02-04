function cowsay(Message,Mood)
% cowsay: cowsay for Matlab
%
% cowsay( [Message] , [Mood] )
%
%  Inputs
% --------
% Message : string, what the cow says (default = 'Hello world !')
% Mood    : string, switch for cow's appearance (default is empty):
%                   -b : borg mode
%                   -d : cow appears dead
%                   -g : greedy mode
%                   -s : stoned cow
%                   -t : tired cow
%                   -y : Young cow
%
%             ~~~~ BONUS MOODS NOT IN ORIGINAL COWSAY ~~~~
%                   -h : happy cow
%                   -o : surprised cow
%
%  Outputs
% ---------
% ;-)
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, June 2018, updated December 2018
%-------------------------------------------------------------------------

if nargin < 1
    Message = 'Hello world !';
end
if nargin < 2
    Mood = '';
end
FLAGmultiline = false;

Message = [' ',Message,' '];

Lmsg = length(Message);

% if message is on multiple lines:
if ~isempty(regexp(Message,'\\n', 'once'))
    Lmsg = max(cellfun(@length,regexp(Message,'\\n','split')));
    Message = vectorize([char(regexp(Message,'\\n','split')'),repmat('\n',length(cellfun(@length,regexp(Message,'\\n','split'))),1)]')';
    Message = Message(1:end-2);
    Message = regexprep(Message,'\\n',' |\\n| ');
    FLAGmultiline = true;
end

SpeechScrollUp = repmat('_',1,Lmsg);
SpeechScrollDown = repmat('-',1,Lmsg);

switch regexprep(lower(Mood),'-','')
    case ''
        CowMood1 = 'oo';
        CowMood2 = ' ';
    case 'b'
        CowMood1 = '==';
        CowMood2 = ' ';
    case 'd'
        CowMood1 = 'xx';
        CowMood2 = 'U';
    case 'g'
        CowMood1 = '$$';
        CowMood2 = ' ';
    case 's'
        CowMood1 = '**';
        CowMood2 = 'U';
    case 't'
        CowMood1 = '--';
        CowMood2 = ' ';
    case 'y'
        CowMood1 = '..';
        CowMood2 = ' ';
    case 'h'
        CowMood1 = '^^';
        CowMood2 = ' ';
    case 'o'
        CowMood1 = 'OO';
        CowMood2 = ' ';
end

if FLAGmultiline
    fprintf([' _',SpeechScrollUp,...
        '_\n| ',Message,...
        ' |\n -',SpeechScrollDown,...
        '-\n        \\   ^__^\n         \\  (',CowMood1,...
        ')\\_______\n            (__)\\       )\\/\\\n             ',CowMood2,...
        '  ||----w |\n                ||     ||\n']);
else
    fprintf([' _',SpeechScrollUp,...
        '_\n< ',Message,...
        ' >\n -',SpeechScrollDown,...
        '-\n        \\   ^__^\n         \\  (',CowMood1,...
        ')\\_______\n            (__)\\       )\\/\\\n             ',CowMood2,...
        '  ||----w |\n                ||     ||\n']);
end

end