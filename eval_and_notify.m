function eval_and_notify(Command)
% Run command and notifies youremail@example.com when execution finishes
% (or whether an error occured (and which kind of error occured))
%
% "Command" is a string with a command such as e.g. :
%   'run /path/to/my/wrapper_script_parallel.m'
%
% If on Windows, requires sendolmail.m and that Outlook is configured.
% If on UNIX system, local mail transfer agent such as sendmail as to be
% configured.
%
% Renaud Marquis @FBMlab, November 2017

if ispc
    
    try
        tic;
        
        eval(Command);
        
        Toc = toc;
        Days = floor(Toc/60/60/24);
        Hours = floor((Toc-Days*60*60*24)/60/60);
        Minutes = floor((Toc-Days*60*60*24-Hours*60*60)/60);
        Seconds = floor((Toc-Days*60*60*24-Hours*60*60-Minutes*60));
        
        sendolmail('youremail@example.com','Command successfully terminated',[Command,' finished in ',num2str(Days),' days, ',num2str(Hours),' hours, ',num2str(Minutes),' minutes and ',num2str(Seconds),' seconds.'])
    catch ME
        if isa(ME,'MException')
            if ~isempty(ME.stack)
                sendolmail('youremail@example.com','Error during command execution',['"',ME.message,'" Line ',num2str(ME.stack(1).line),' in ',ME.stack(1).file])
                warning(['"',ME.message,'" Line ',num2str(ME.stack(1).line),' in ',ME.stack(1).file])
            else
                sendolmail('youremail@example.com','Error during command execution',ME.message)
                warning(ME.message)
            end
        end
    end
    
elseif isunix
    
    try
        tic;
        
        eval(Command);
        
        Toc = toc;
        Days = floor(Toc/60/60/24);
        Hours = floor((Toc-Days*60*60*24)/60/60);
        Minutes = floor((Toc-Days*60*60*24-Hours*60*60)/60);
        Seconds = floor((Toc-Days*60*60*24-Hours*60*60-Minutes*60));
        
        system('echo "',[Command,' finished in ',num2str(Days),' days, ',num2str(Hours),' hours, ',num2str(Minutes),' minutes and ',num2str(Seconds),' seconds.'],'" | mail -s "Command successfully terminated" "youremail@example.com"')
    catch ME
        if isa(ME,'MException')
            if ~isempty(ME.stack)
                system('echo "',['"',ME.message,'" Line ',num2str(ME.stack(1).line),' in ',ME.stack(1).file],'" | mail -s "Error during command execution" "youremail@example.com"')
                warning(['"',ME.message,'" Line ',num2str(ME.stack.line),' in ',ME.stack.file])
            else
                system('echo "',ME.message,'" | mail -s "Error during command execution" "youremail@example.com"')
                warning(ME.message)
            end
        end
    end
    
else
    
    error('Don''t know what kind of operating system this is...')
    
end
