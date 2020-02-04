function [FilePaths,DirPaths] = my_super_listfiles(Directory,RegexFilterIn,Recursive,varargin)
% my_super_listfiles: this is just an alias for spm_select.m that avoids
% having to look at syntax again, and which allows more powerful regular
% expressions that my_listfiles and my_recursive_listfiles.
%
% [ FilePaths, DirPaths ] = my_super_listfiles( Directory, RegexFilterIn )
% [ FilePaths, DirPaths ] = my_super_listfiles( Directory, RegexFilterIn, Recursive )
% [ FilePaths, DirPaths ] = my_super_listfiles( Directory, RegexFilterIn, Recursive, RegexFilterOut, RegexFilterOutLevel )
%
%  Inputs
% --------
% Directory:                 char, directory in which to look for
%                            files / directories
% RegexFilterIn:             char, regular expression to match files, which
%                            supports things such as e.g. :
%      '^customprefix.*(basename1|basename2).*customsuffix.*\.customextension$'
% Recursive [optional]:      logical, whether to look in child directories
%                            of Directory (1, default) or not (0)
%
% -------------- FOR POWER USERS --------------
% RegexFilterOut [optional]: char, regular expression to filter out (or filter
%                            in) files, which also supports regex like RegexFilterIn
% RegexFilterOutLevel [optional] : char, indicating the level at which
%                            RegexFilterOut should operate, among
%                            possible options 'path', 'basename', 
%                            'ext', 'filename', 'cpath' or 'fpath'.
%                            See spm_file.m for more details.
% RegexFilterInOrOut [optional]: char, indicating whether RegexFilterOut
%                            should exclude ('out') or keep ('in') files /
%                            paths / extensions, ... etc.
%
% RegexFilterOut and RegexFilterOutLevel should be introduced as triplets,
% and can be multiple triplets of RegexFilterOut, RegexFilterOutLevel, and 
% RegexFilterInOrOut, where each RegexFilterOut applies at different levels
% ('ext', 'basename', 'path'), e.g.:
%           [ FilePaths, DirPaths ] = ...
%               my_super_listfiles( '/my/root/path', '^hdEEG.*', 1, ...
%               'ext', 'sef', 'in',...
%               'path', 'Epochs', 'in',...
%               'basename', '(interp|filtered)', 'out',...
%               'path', 'More', 'out')
%
%  Outputs
% ---------
% FilePaths:   char array, absolute path(s) to file(s)
% DirPaths:    char array, absolute path(s) to directory(ies) containing
%              FilePaths
%
% See also MY_LISTFILES, MY_RECURSIVE_LISTFILES, SPM_FILE, SPM_SELECT
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, October 2019
%-------------------------------------------------------------------------

if nargin < 3
    Recursive = true;
end

if Recursive
    FilePaths = spm_select('FPListRec',Directory,RegexFilterIn);
else
    FilePaths = spm_select('FPList',Directory,RegexFilterIn);
end

if nargin > 3
    if mod(length(varargin),3)~=0
        error('RegexFilterOut and RegexFilterOutLevel must be provided as triplets')
    else
        FilterOutMatrix = reshape(varargin,3,length(varargin)/3);
    end
    
    Options = {'path','basename','ext','filename','cpath','fpath'};
    
    for f = 1:size(FilterOutMatrix,2)
        RegexFilterOut = FilterOutMatrix{2,f};
        FilterAt = FilterOutMatrix{1,f};
        FilterInOrOut = FilterOutMatrix{3,f};
        
        if ismember(FilterAt,Options)
            switch lower(FilterInOrOut)
                case 'out'
                    FilePaths = filter_out(FilePaths,RegexFilterOut,FilterAt);
                case 'in'
                    FilePaths = filter_in(FilePaths,RegexFilterOut,FilterAt);
                otherwise
                    warning('Filter level "%s" is not a recognized options, ignoring...',FilterInOrOut)
            end
        else
            warning('Filter level "%s" is not a recognized options, ignoring...',FilterAt)
        end
    end
    
    DirPaths = spm_file(FilePaths,'path');
else
    DirPaths = spm_file(FilePaths,'path');
end

end

function FilePaths = filter_out(FilePaths,RegexFilterOut,FilterAt)

ExcludeThese = cellfun(@isempty,regexp(cellstr(spm_file(FilePaths,FilterAt)),RegexFilterOut));
FilePaths = cellstr(FilePaths);
FilePaths = FilePaths(ExcludeThese);
FilePaths = char(FilePaths);

end

function FilePaths = filter_in(FilePaths,RegexFilterOut,FilterAt)

IncludeThese = ~cellfun(@isempty,regexp(cellstr(spm_file(FilePaths,FilterAt)),RegexFilterOut));
FilePaths = cellstr(FilePaths);
FilePaths = FilePaths(IncludeThese);
FilePaths = char(FilePaths);

end
