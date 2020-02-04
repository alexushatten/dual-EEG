function mkr = dual_read_mkr(txtfile,DINtype)

% NB: need to add .txt to the file in finder
% NB: may need to edit manually marker file e.g. 256(2)
% GOAL      read markers from text files for dual EEG recordings
% DINtype: 'DINin' or 'DINafter' (whether the number to extract is "128" or
%          "12" in e.g. "DIN12:128")
% INPUT     text file with time stamps (column 1 & 2) and marker labels (column 3)
% OUTPUT    mkr is a two column matrix: first column is trigger onset in datapoints, second
%           column is marker label

% LOG:
%   - developped March 2017, Maxime Baud
%   - July 2017 changed regexp function to be able to read all marker files
%   - refacto, Renaud Marquis, August 2017, last updated July 2019

mkr = do_read_mkr(txtfile,'num');

if isempty(mkr)
    warning('NO NUMERIC MARKERS FOUND IN %s',txtfile)
    warning('Trying to trim labels with pattern "DIN"...')
    try % #RM@FBMlab, 2019-07-01: to handle cases without numeric markers more explicitly
        mkr = do_read_mkr(txtfile,DINtype);
    catch % #RM@FBMlab, 2019-07-01 
        error('NO NUMERIC MARKERS FOUND IN %s',txtfile) % #RM@FBMlab, 2019-07-01 
    end % #RM@FBMlab, 2019-07-01 
    if isempty(mkr)
        error('NO NUMERIC MARKERS FOUND IN %s',txtfile)
    end
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

function mkr = do_read_mkr(txtfile,opt)
% txtfile: path to .mrk file
% opt: 'num' for numeric markers expected ; 'DIN' otherwise

fid = fopen(txtfile);
i=0;
while (feof(fid)==0)    % while  not at the end of the file
    line = fgetl(fid); %take line by line
    pieces = regexp(line,'\w+','match');
    for p=1:numel(pieces)
        pcs(1,p)  = str2double(pieces(p));
    end
    if ~isnan(pcs(1)) % skip header
        if strcmpi(opt,'num')
            if numel(pcs)>2 && ~isnan(pcs(3)) % label must be a number for it to be a trigger
                i=i+1;
                mkr(i,1) = str2double(pieces{1}); % trigger onset
                mkr(i,2) = str2double(pieces{3}); % trigger label
                %         pieces = regexp(line, '\t', 'split');  %match  \t (tab) and split into an array of cells
                %         label = regexp(pieces{3},'"','split');
                %         mkr(i,2) = str2double(label{2
            end
        elseif strcmpi(opt,'DINafter')
            if numel(pcs)>3 && strcmp(pieces{3}(1:3),'DIN') % based on  ALB TPG1-TPG2
                i=i+1;
                mkr(i,1) = str2double(pieces{1}); % trigger onset
                mkr(i,2) = str2double(pieces{4}); % trigger label NOT WHAT IS RIGHT AFTER "DIN" BUT WHAT IS AFTER ":" (THUS THE 4TH ELEMENT)
                %         pieces = regexp(line, '\t', 'split');  %match  \t (tab) and split into an array of cells
                %         label = regexp(pieces{3},'"','split');
                %         mkr(i,2) = str2double(label{2
            end
        elseif strcmpi(opt,'DINin')
            if numel(pcs)>3 && strcmp(pieces{3}(1),'D') % based on CC rest part 1 subpart 2
                i=i+1;
                mkr(i,1) = str2double(pieces{1}); % trigger onset
                NumIdx = regexp(pieces{3},'[0-9]');
                mkr(i,2) = str2double(pieces{3}(NumIdx));
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