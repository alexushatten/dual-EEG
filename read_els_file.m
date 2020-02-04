function [x,y,z,name,ClusterName,Nclus,FullName,Els] = read_els_file(Filename)
% READ_ELS_FILE: reads Cartool .els file
%
% Usage:
%-------------------------------------------------------------------------
%
% [x,y,z,name,ClusterName,Nclus,FullName,Els] = read_els_file(Filename)
%
% Outputs:
%-------------------------------------------------------------------------
%
% x, y, z coordinates
%
% name: name of electrode / contact
%
% ClusterName: name of cluster (group of electrodes / contacts)
%
% Nclus: number of clusters
%
% Els: the whole text file as a data array
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, February 2018, last update November 2018
%-------------------------------------------------------------------------

formatSpec = '%s%[^\n\r]';
fileID = fopen(Filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '',  'ReturnOnError', false);
dataArray{1} = strtrim(dataArray{1});
fclose(fileID);
Els = dataArray{:, 1};

if ~strcmp(Els(1),'ES01')
    warning('Maybe not an .els file! (magic number mismatch)')
end

PresNelec = regexp(Els(2,:),'\d','match');
if ~isempty(PresNelec{1})
    PresNelec = str2double(cell2mat(PresNelec{1}));
end
PresNclus = regexp(Els(3,:),'\d','match');
if ~isempty(PresNclus{1})
    PresNclus = str2double(cell2mat(PresNclus{1}));
end
name = {};
x = [];
y = [];
z = [];
Nclus = 0;
NextMightBeNcon = 0;
Count = 0;
WarnBefore = false;
for l = 4:size(Els,1)
    if NextMightBeNcon == 1
        Ncon = regexp(Els(l,:),'\d','match');
        Ncon = str2double(cell2mat(Ncon{1}));
        NextMightBeNcon = 0;
    else
        Parsed = regexp(Els(l,:),'\s|\t','split');
        if size(Parsed{1},2)>2
            % it is an electrode, with coordinates and name
%             if size(Parsed{1},2)==1
%                 Dim{Count+1}=str2double(Parsed{1}{1});
% %             elseif size(Parsed{1},2)~=4
%             elseif size(Parsed{1},2)>1
            if size(Parsed{1},2)>1
                if size(Parsed{1},2)==3
                    CountCon = CountCon+1;
                    Count = Count+1;
                    x(Count) = str2double(Parsed{1}{1});
                    y(Count) = str2double(Parsed{1}{2});
                    z(Count) = str2double(Parsed{1}{3});
                    name{Count} = num2str(CountCon);
                    FullName{Count} = [ClusterName{Nclus},num2str(CountCon)];
                else
                    % ==============================================
                    %  try this as a last resort:
                    Parsed = regexp(Els(l,:),'\t|\s','split');
                    
                    if size(Parsed{1},2)>2
                        % it is an electrode, with coordinates and name
                        if size(Parsed{1},2)~=4
                            if ~WarnBefore
                                warning('Do all electrodes/contacts have a proper name?')
                                WarnBefore = true;
                            end
                            if size(Parsed{1},2)==3
                                CountCon = CountCon+1;
                                Count = Count+1;
                                x(Count) = str2double(Parsed{1}{1});
                                y(Count) = str2double(Parsed{1}{2});
                                z(Count) = str2double(Parsed{1}{3});
                                name{Count} = num2str(CountCon);
                                FullName{Count} = [ClusterName{Nclus},num2str(CountCon)];
                            else
                                error('Can''t parse .els file...')
                            end
                        else
                            Count = Count+1;
                            x(Count) = str2double(Parsed{1}{1});
                            y(Count) = str2double(Parsed{1}{2});
                            z(Count) = str2double(Parsed{1}{3});
                            name{Count} = Parsed{1}{4};
                            FullName{Count} = [ClusterName{Nclus},Parsed{1}{4}];
                        end
                    else
                        Test = regexp(Els(l,:),'\D','match');
                        if ~isempty(Test{1})
                            Nclus = Nclus+1;
                            Temp = regexp(Els(l,:),'\D','match');
                            ClusterName{Nclus} = char(Temp{1})';
                            CountCon = 0; % reset it in case contacts do not have specific name
                            NextMightBeNcon = 1;
                        end
                    end
                end
            end
        else
            Test = regexp(Els(l,:),'\D','match');
            if ~isempty(Test{1})
                Nclus = Nclus+1;
                Temp = regexp(Els(l,:),'\D','match');
                ClusterName{Nclus} = char(Temp{1})';
                CountCon = 0; % reset it in case contacts do not have specific name
                NextMightBeNcon = 1;
            end
        end
    end
end

if PresNclus~=Nclus
    warning('Header of .els possibly wrong, or some groups or electrodes/contacts are missing')
end
if size(x,2)~=PresNelec;
    warning('Header of .els possibly wrong, or some electrodes contacts are missing')
end

x = x';
y = y';
z = z';
name = name';

end
