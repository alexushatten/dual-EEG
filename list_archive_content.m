function [Out, Unfiltered, OutAll] = list_archive_content( Path, Filter, DateRange )
% list_archive_content: lists archive content using 7zip
%
%  Inputs
% --------
% Path   : string, path to archive file, e.g. '/path/to/archive.zip'
%
% Filter : string, regex for filtering whose wildcard characters will
%   automatically be escaped with regexptranslate('escape',Filter). If
%   DateRange is set to true, set it to YYYY-MM-DD format and use semicolon
%   to separate dates, like this: YYYY-MM-DD:YYYY-MM-DD
% 
% DateRange : flag set to false by default, set it to true to search for
%   files within a certain date range.
%
%-------------------------------------------------------------------------
% https://www.7-zip.org/
%
% /!\ Define path to 7-Zip in "PathTo7z"
%
% Nota bene:
%  On UNIX systems, you can specify the path to 7z executable, but you can
%  also simply use "zipinfo", "unzip -l", or "lsar" in combination with grep:
%   $ zipinfo 'Path' | grep 'Filter'
%
%  The Filter applies to filename, path, extension, date and file
%  attribute, so you can choose to display only files using:
%   >> list_archive_content(Path,'....A')
%
%  However, if you want to search for files within a certain date range,
%  please use the flag DateRange and specify Filter as YYYY-MM-DD:YYYY-MM-DD.
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, June 2018
%-------------------------------------------------------------------------

% Locate 7-zip:
PathTo7z = 'C:\Program Files\7-Zip\7z.exe';
% PathTo7z = '';
P27z = '';

if nargin < 3
    DateRange = false;
end
if DateRange && (isempty(regexp(Filter,':','once')) || length(regexp(Filter,':'))>1)
    error('Improper format of Filter with DateRange enabled, check help')
end
if nargin < 2
    Filter = '';
end

% Locate 7-zip:
if ispc
    if isempty(PathTo7z)
        [~,P27z] = system('where /r "C:\Program Files" 7z');
    end
    if isempty(strtrim(P27z)) && isempty(PathTo7z)
        [~,P27z] = system('where /r "C:\Program Files (x86)" 7z');
        PathTo7z = strtrim(P27z);
    end
else
    [~,PathTo7z] = system('which 7z');
end
if isempty(PathTo7z) || strcmp(PathTo7z,'7z not found')
    error('Path to 7-Zip not found')
end

% Command to list all content from archive file
[~,Unfiltered] = system(['"',PathTo7z,'" l ','"',Path,'"']);

% Split output according to newline (\n)
if exist('splitlines','file')~=2
    error('Matlab version < R2016b, please get alternative function "splitLines.m" for converting 7-zip output to array')
else
    [~,Check] = fileparts(which('splitlines'));
end
if strcmp(Check,lower(Check))
    OutA = splitlines(Unfiltered)';
else
    eval(['OutA = ',Check,'(Unfiltered);']);
    OutA = OutA'; %#ok<NODEF>
end

% Filter
if isempty(Filter)
    KeepEm = true(size(OutA));
else
    Table = find(~cellfun(@isempty,regexp(OutA,'-------------------')));
    OutA = OutA(Table(1)+1:Table(2)-1);
    if DateRange
        Dates = datenum(regexp(Filter,':','split'),'YYYY-MM-DD');
        AllDates = cell(Table(2)-Table(1)-1,1);
        for d = 1:size(OutA,1)
            TempD = regexp(OutA{d},' ','split');
            AllDates{d} = TempD{1};
        end
        AllDatesNum = datenum(AllDates,'YYYY-MM-DD');
        KeepEm = ((AllDatesNum >= min(Dates)) ++ (AllDatesNum <= max(Dates)))>1;
    else
        
        KeepEm = ~cellfun(@isempty,regexp(OutA,regexptranslate('escape',Filter)));
    end
end

Out = cell(size(OutA(KeepEm))+[2 0]);
Out{1} = '   Date      Time    Attr         Size   Compressed  Name ';
Out{2} = '------------------- ----- ------------ ------------  ------------------------ ';
Out(3:end) = OutA(KeepEm);
Out = char(Out);

OutAll = OutA;

end

