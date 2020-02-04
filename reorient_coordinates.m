function CoordsA = reorient_coordinates( Coords, A )
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

% CoordsA = [Coords(:,1) Coords(:,2) Coords(:,3) ones(size(Coords,1),1)]*(inv(A))';
CoordsA = A*[Coords(:,1) Coords(:,2) Coords(:,3) ones(size(Coords,1),1)]';
CoordsA = CoordsA';
CoordsA(:,4) = [];

end

