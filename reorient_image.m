function reorient_image(path2image,A)
%reorient_coordinates: apply affine transformation matrix to coordinates
%
% Usage:
%-------------------------------------------------------------------------
% 
% CoordsA = reorient_coordinates( Coords, A )
%
% Inputs:
%-------------------------------------------------------------------------
%
% Coords: n x 3 coordinates matrix
%
% A: affine transformation matrix given e.g. by spm_matrix
%
% Outputs:
%-------------------------------------------------------------------------
%
% CoordsA: updated coordinates
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, February 2018
%-------------------------------------------------------------------------


matlabbatch{1}.spm.util.reorient.srcfiles = {[path2image,',1']};
matlabbatch{1}.spm.util.reorient.transform.transM = A;
matlabbatch{1}.spm.util.reorient.prefix = '';

spm_jobman('run',matlabbatch)

end
