function [mkr,IsNumeric] = dual_read_mkr_NEW(txtfile)

% NB: need to add .txt to the file in finder
% NB: may need to edit manually marker file e.g. 256(2)
% GOAL      read markers from text files for dual EEG recordings
% INPUT     text file with time stamps (column 1 & 2) and marker labels (column 3)
% OUTPUT    mkr is a two column matrix: first column is trigger onset in datapoints, second
%           column is marker label

% LOG:
%   - developped March 2017, Maxime Baud
%   - July 2017 changed regexp function to be able to read all marker files
%   - refacto, Renaud Marquis, August 2017, last updated July 2019

[mkr,IsNumeric] = do_read_mkr(txtfile,'num');

if isempty(mkr)
    warning('Trying to trim labels with pattern "DIN"...')
    try % #RM@FBMlab, 2019-07-01: to handle cases without numeric markers more explicitly
        [mkr,IsNumeric] = do_read_mkr(txtfile,'DIN');
    catch % #RM@FBMlab, 2019-07-01
        warning('Trying to trim labels with pattern "EprimeNetStation_DIN"...')
        try
            [mkr,IsNumeric] = do_read_mkr(txtfile,'EprimeNetStation_DIN');
        catch
            warning('Trying to trim labels with pattern "EprimeMicromed_"...')
            try
                [mkr,IsNumeric] = do_read_mkr(txtfile,'EprimeMicromed_');
            catch
                warning('Trying to trim labels with pattern "EprimeNetStation_"...')
                try
                    [mkr,IsNumeric] = do_read_mkr(txtfile,'EprimeNetStation_');
                catch
                    error('NO NUMERIC MARKERS FOUND IN %s',txtfile) % #RM@FBMlab, 2019-07-01
                end
            end
        end
    end % #RM@FBMlab, 2019-07-01
    if isempty(mkr)
        warning('Trying to trim labels with pattern "EprimeNetStation_DIN"...')
        %         try
        [mkr,IsNumeric] = do_read_mkr(txtfile,'EprimeNetStation_DIN');
        %         catch
        if isempty(mkr)
            warning('Trying to trim labels with pattern "EprimeMicromed_"...')
            %             try
            [mkr,IsNumeric] = do_read_mkr(txtfile,'EprimeMicromed_');
            %             catch
            if isempty(mkr)
                warning('Trying to trim labels with pattern "EprimeNetStation_"...')
                %                 try
                if isempty(mkr)
                    [mkr,IsNumeric] = do_read_mkr(txtfile,'EprimeNetStation_');
                    %                 catch
                    if isempty(mkr)
                        error('NO NUMERIC MARKERS FOUND IN %s',txtfile) % #RM@FBMlab, 2019-07-01
                    end
                end
            end
        end
    end
%                 if isempty(mkr)
%                     error('NO NUMERIC MARKERS FOUND IN %s',txtfile) % #RM@FBMlab, 2019-07-01
%             end
%         end
%     end
end

% Inform on markers
if mkr(1,2)~=1
    warning('Missing markers at the beginning')
else warning('No missing marker at the beginning')
end
if any(unique(diff(mkr(:,2)))~=1)
    warning('Missing markers in the middle')
else warning('No missing markers in the middle')
end

end

function [mkr,IsNumeric] = do_read_mkr(txtfile,opt)
% txtfile: path to .mrk file
% opt: 'num' for numeric markers expected ; 'DIN' otherwise

fid = fopen(txtfile);
i=0;
Count = 0;
while (feof(fid)==0)    % while  not at the end of the file
    line = fgetl(fid); %take line by line
    pieces = regexp(line,'\w+','match');
    for p=1:numel(pieces)
        pcs(1,p)  = str2double(pieces(p));
    end
    if ~isnan(pcs(1)) % skip header
        Count = Count+1;
        if strcmpi(opt,'num')
            if numel(pcs)>2 && ~isnan(pcs(3)) % label must be a number for it to be a trigger
                i=i+1;
                mkr(i,1) = str2double(pieces{1}); % trigger onset
                mkr(i,2) = str2double(pieces{3}); % trigger label
                %         pieces = regexp(line, '\t', 'split');  %match  \t (tab) and split into an array of cells
                %         label = regexp(pieces{3},'"','split');
                %         mkr(i,2) = str2double(label{2
                IsNumeric(Count) = true;
            else
                IsNumeric(Count) = false;
            end
        elseif strcmpi(opt,'DIN')
            if numel(pcs)>3 && strcmp(pieces{3}(1:3),'DIN') % based on patient1
                i=i+1;
                mkr(i,1) = str2double(pieces{1}); % trigger onset
                mkr(i,2) = str2double(pieces{4}); % trigger label NOT WHAT IS RIGHT AFTER "DIN" BUT WHAT IS AFTER ":" (THUS THE 4TH ELEMENT)
                %         pieces = regexp(line, '\t', 'split');  %match  \t (tab) and split into an array of cells
                %         label = regexp(pieces{3},'"','split');
                %         mkr(i,2) = str2double(label{2
                IsNumeric(Count) = true;
            else
                IsNumeric(Count) = false;
            end
        elseif strcmpi(opt,'EprimeNetStation_')
            if ~isempty(regexp(pieces{3},'EprimeNetStation_','once'))
%             if numel(pcs)==3 && strcmp(pieces{3}(1:17),'EprimeNetStation_') % based on patient1
                i=i+1;
                mkr(i,1) = str2double(pieces{1}); % trigger onset
                mkr(i,2) = str2double(regexprep(pieces{3},'EprimeNetStation_','')); % trigger label NOT WHAT IS RIGHT AFTER "DIN" BUT WHAT IS AFTER ":" (THUS THE 4TH ELEMENT)
                %         pieces = regexp(line, '\t', 'split');  %match  \t (tab) and split into an array of cells
                %         label = regexp(pieces{3},'"','split');
                %         mkr(i,2) = str2double(label{2
                IsNumeric(Count) = true;
            else
                IsNumeric(Count) = false;
            end
        elseif strcmpi(opt,'EprimeNetStation_DIN')
%             if numel(pcs)>3 && strcmp(pieces{3}(1:20),'EprimeNetStation_DIN') % based on patient1
            if ~isempty(regexp(pieces{3},'EprimeNetStation_DIN','once'))
                i=i+1;
                mkr(i,1) = str2double(pieces{1}); % trigger onset
                mkr(i,2) = str2double(regexprep(pieces{3},'EprimeNetStation_DIN','')); % trigger label NOT WHAT IS RIGHT AFTER "DIN" BUT WHAT IS AFTER ":" (THUS THE 4TH ELEMENT)
                %         pieces = regexp(line, '\t', 'split');  %match  \t (tab) and split into an array of cells
                %         label = regexp(pieces{3},'"','split');
                %         mkr(i,2) = str2double(label{2
                IsNumeric(Count) = true;
            else
                IsNumeric(Count) = false;
            end
        elseif strcmpi(opt,'EprimeMicromed_')
            if ~isempty(regexp(pieces{3},'EprimeMicromed_','once'))
                i=i+1;
                mkr(i,1) = str2double(pieces{1}); % trigger onset
                mkr(i,2) = str2double(regexprep(pieces{3},'EprimeMicromed_','')); % trigger label NOT WHAT IS RIGHT AFTER "DIN" BUT WHAT IS AFTER ":" (THUS THE 4TH ELEMENT)
                %         pieces = regexp(line, '\t', 'split');  %match  \t (tab) and split into an array of cells
                %         label = regexp(pieces{3},'"','split');
                %         mkr(i,2) = str2double(label{2
                IsNumeric(Count) = true;
            else
                IsNumeric(Count) = false;
            end
        else
            error('CHECK CODE, unknown option %s for subfunction do_read_mkr',opt)
        end
    end
end
fclose(fid);

if exist('mkr')~=1 % (and not just ~... as could be > 1)
    mkr = [];
end

end