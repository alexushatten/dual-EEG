function matlab2mgrid_RM(mgridFname,NewCoords,Subtract1)
% MATLAB2MGRID_RM:  replace coordinates in .mgrid files with updated ones
%
%
% Usage:
%-------------------------------------------------------------------------
% matlab2mgrid_RM(mgridFname,NewCoords,Subtract1)
%
%
% Inputs:
%-------------------------------------------------------------------------
% mgridFname    : string, path to .mgrid file 
%
% NewCoords     : n x 3 array of 3D coordinates (double)
%
% Subtract1     : BioImageSuite indexes the first voxel as [0 0 0] but
%                 mgrid2matlab.m adds 1 to coordinates. If this has not been
%                 corrected elsewhere, default [1] is to subtract 1 to
%                 NewCoords before inserting them into new .mgrid file.
%
%
% Outputs:
%-------------------------------------------------------------------------
% New .mgrid file (mgridFname with "_align_equ" appended) with NewCoords
% replacing actual coordinates of electrodes
%
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, January 2018
%-------------------------------------------------------------------------

if nargin < 3
    Subtract1 = 1;
end

if Subtract1==1
    NewCoords = NewCoords - 1; % BioImageSuite indexes the first voxel as [0 0 0]
end

[p,n,e] = fileparts(mgridFname);

NewMgridFname = strcat(p,filesep,[n,'_align_equ',e]);

% Read .mgrid file
fid = fopen(mgridFname,'r');
i = 1;
tline = fgetl(fid);
A{i} = tline;
Cnt = 0;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    
    % If previous line declared #Position, insert new coordinates instead
    % of actual ones
    if strcmpi(A{i-1}, '#Position')
        Cnt = Cnt + 1;
        A{i} = sprintf('% .4f',NewCoords(Cnt,:));
    else
        A{i} = tline;
    end
    
end
fclose(fid);

if Cnt ~= size(NewCoords,1), error('Size of coordinates matrices don''t match, exiting...'), end;

% Write new .mgrid file with updated coordinates
fid = fopen(NewMgridFname,'w');
for i = 1:numel(A)
    if A{i+1} == -1
        fprintf(fid,'%s', A{i});
        break
    else
        fprintf(fid,'%s\n', A{i});
    end
end
fclose(fid);

end

