function write_els_file(FileName,xyz,name,fullnames,electype)
%write_els_file: does NOT work yet with grids, only depth and strips
%
% Usage:
%-------------------------------------------------------------------------
% 
% write_els_file(FileName,xyz,name)
%
% Inputs:
%-------------------------------------------------------------------------
%
% FileName: filename, with extension (written in pwd if not full file path)
%
% xyz: [n x 3] coordinates matrix
%
% name: [n x 1] electrode labels
%
% fullnames: flag, whether to insert contacts names (1) or just number within
%           electrodes (0)
%
% electype: [n x 1] whether each electrode is part of a depth/strip (1), grid (2)
%           or head cap (3), see Cartool reference manual for more details
%
%-------------------------------------------------------------------------
% NB: channels have to be grouped by electrode type (for a depth electrode
% named "ABC", contacts within "ABC" have to follow each other)!
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, February 2018
%-------------------------------------------------------------------------

if nargin<5
    electype = ones(length(name),1);
end

fileID = fopen(FileName,'w');
fprintf(fileID,'%3s\r\n','ES01');

[~,ElecMatch,~,~] = bipolar_montage(name);

clear ChNames
for c = 1:length(name)
    ChNames{c} = char(regexpi(name{c},'[a-z]+','match')); %#ok<AGROW>
end
ChNames = ChNames';
ElecNames = unique(ChNames,'stable');

NcontactsPerElec = diff([0;find(diff(ElecMatch)~=0);size(xyz,1)]);

if size(ElecNames,1)~=size(NcontactsPerElec,1),error('Some electrodes seem to have the same label?!?^'),end

Count = 0;
% how many contacts in total
fprintf(fileID,'%d\r\n',sum(NcontactsPerElec));
% how many electrodes (groups of contacts)
fprintf(fileID,'%d\r\n',size(NcontactsPerElec,1));
for elec = 1:size(NcontactsPerElec,1)
    
    % write electrode name
    fprintf(fileID,'%s\r\n',upper(ElecNames{elec}));
    
    % write number of contacts
    fprintf(fileID,'%d\r\n',NcontactsPerElec(elec));
    
    if electype(Count+1) == 0 || electype(Count+1) == 1 % for lines and points
        % write "1"
        fprintf(fileID,'%d\r\n',1);
    elseif electype(Count+1) == 2 % for 2D setup (like grids)
        fprintf(fileID,'%d\r\n',2);
    elseif electype(Count+1) == 3 % for 3D setup (like head cap)
        fprintf(fileID,'%d\r\n',3);
    else
        error('Unrecognized electrode type!')
    end
    
    for cont = 1:NcontactsPerElec(elec)
        Count = Count+1;
        
        if fullnames==1
            formatSpec = '\t%.7f %.7f %.7f %s\r\n';
            fprintf(fileID,formatSpec,xyz(Count,1),xyz(Count,2),xyz(Count,3),name{Count});
        elseif fullnames==0
            formatSpec = '\t%.7f %.7f %.7f %s\r\n';
%             fprintf(fileID,formatSpec,xyz(Count,1),xyz(Count,2),xyz(Count,3),cont);
%             fprintf(fileID,formatSpec,xyz(Count,1),xyz(Count,2),xyz(Count,3),NcontactsPerElec(elec)-(cont-1)); % should be done in reverse order (as done in BioImageSuite and as recommended by iELVis developpers)
            fprintf(fileID,formatSpec,xyz(Count,1),xyz(Count,2),xyz(Count,3),char(regexpi(name{Count},'[0-9]+','match')));
        end
    end
    
end
fclose(fileID);

end