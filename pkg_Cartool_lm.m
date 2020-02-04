function pkg_Cartool_lm( LMfile, TargetPath, TargetTarget )
% pkg_Cartool_lm: package a Cartool .lm (Link Many) file, i.e. gets all
% files referenced and put them in a target directory
%
% pkg_Cartool_lm( LMfile, TargetPath )
%
%  Inputs
% --------
% LMfile:       path to .lm file
% TargetPath:   path where to rewrite all files referenced in original .lm
%               file
% TargetTarget: [optional] path that should be used within new .lm file if
%              it is planned to be used elsewhere than locally
%
%  Outputs
% ---------
% Directory @ TargetPath containing all files
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, March 2019
%-------------------------------------------------------------------------

formatSpec = '%s%[^\n\r]';
fileID = fopen(LMfile,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '',  'ReturnOnError', false);
dataArray{1} = strtrim(dataArray{1});
fclose(fileID);
LMcontent = dataArray{:, 1};

if ~(exist(TargetPath,'file')==7)
    mkdir(TargetPath);
end

for f = 1:length(LMcontent)
    if ~(exist(LMcontent{f},'file')==2)
        warning('File %s could not be found, please select it',LMcontent{f})
        LMcontent{f} = spm_select;
    end
    
    copyfile(LMcontent{f},TargetPath)
end

% self-copy & change path

[~,LMfilename,LMext] = fileparts(LMfile);
TargetLM = fullfile(TargetPath,[LMfilename,LMext]);

if nargin > 2
    TargetPath = TargetTarget;
end

if ~isempty(regexp(TargetPath,'\', 'once')) && ~isempty(regexp(TargetPath,'/', 'once'))
    error('Target path contains both "/" and "\", does not know how to parse it!')
elseif isempty(regexp(TargetPath,'/', 'once'))
    TargetIsUNIX = false;
elseif isempty(regexp(TargetPath,'\', 'once'))
    warning('Cartool on UNIX system?! Are you kidding? Okay, processing anyway...')
    TargetIsUNIX = true;
else
    error('Cannot parse target path "%s"',TargetPath)
end

NewLMcontent = cell(size(LMcontent));

for l = 1:length(LMcontent)
    [~,asdf1,asdf2] = fileparts(LMcontent{l});
    FileNameExt = [asdf1,asdf2];
    if TargetIsUNIX
        NewLMcontent{l} = regexprep(fullfile(TargetPath,FileNameExt),'\','/');
    else
        NewLMcontent{l} = regexprep(fullfile(TargetPath,FileNameExt),'/','\');
    end
end

fileID = fopen(TargetLM,'w');
for l = 1:length(NewLMcontent)
    fprintf(fileID,'%s\r\n',char(NewLMcontent{l}));
end
fclose(fileID);

end

