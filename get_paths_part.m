function FolderLevel = get_paths_part(Paths,Level)
% get_paths_part: parses paths and get i-th folder level from root
%
% FolderLevel = get_paths_part(Paths)
%
%  Inputs
% --------
% Paths:       char array of file paths
% Level:       integer, folder level from root
%
%  Outputs
% ---------
% FolderLevel: cell of char array of file paths with everything stripped
%              except folder "Level"
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, November 2019
%-------------------------------------------------------------------------

if ispc
    Parsed = regexp(Paths,repmat(filesep,1,2),'split'); % self-escape needed on Windows
else
    Parsed = regexp(Paths,filesep,'split');
end

FolderLevel = cell(size(Parsed));
for fp = 1:size(Parsed,1)
    FolderLevel(fp) = Parsed{fp}(Level);
end

end

