function mkr = dual_read_mkr(txtfile)

% NB: need to add .txt to the file in finder
% GOAL      read markers from text files for dual EEG recordings
% INPUT     text file with datapoints (column 1) and marker label (column 3)
% OUTPUT    mkr is a two column matrix: first column is datapoints, second
%           column is marker label


fid = fopen(txtfile);
i=0;
while (feof(fid)==0)    % while  not at the end of the file
    line = fgetl(fid); %take line by line
    if line(1)==' '; % skip header
        i=i+1;
        pieces = regexp(line, '\t', 'split');  %match  \t (tab) and split into an array of cells
        mkr(i,1) = str2double(pieces{1}); % fill in matrix
        label = regexp(pieces{3},'"','split');
        mkr(i,2) = str2double(label{2});
    end
end
fclose(fid); %#RM, 04.07.17 11:32:56
