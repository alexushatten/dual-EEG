function xyzhat_eq = AlignEqualizeContactDepthElec(xyz)
% AlignEqualizeContactDepthElec:  align 3D coordinates of contacts on depth
%                       electrodes and make distance between contacts equal
%
% Usage:
%-------------------------------------------------------------------------
% xyzhat_eq = AlignEqualizeContactDepthElec(xyz)
%
% Inputs:
%-------------------------------------------------------------------------
% xyz :         N x 3 matrix of coordinates
%
% Outputs:
%-------------------------------------------------------------------------
% xyzhat_eq :   N x 3 matrix of coordinates after alignment and interpolation
%
% Nota bene:
%-------------------------------------------------------------------------
% - Assumes first and last contacts are well positionned
% - Uses SVD to find principal vector
%
% Acknowledgements:
%-------------------------------------------------------------------------
% - This function makes use of David Legland's geom3d
% http://mathworks.com/matlabcentral/fileexchange/24484-geom3d
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, January 2018
%-------------------------------------------------------------------------

% % find principal vector
% r0=mean(xyz);
% xyz=bsxfun(@minus,xyz,r0);
% [~,~,V]=svd(xyz,0);
% 
% % make vector and project original coordinates on it
% R3D = (kron(1:2,V(:,1))+repmat(r0',1,2))';
% xyzhat = projPointOnLine3d(xyz,[R3D(1,:),R3D(2,:)-R3D(1,:)]);

% find principal vector and project points back on line:
L3D = fitLine3d(xyz);
xyzhat = projPointOnLine3d(xyz,L3D);

% equalize distance between contacts (assumes first and last contacts are
% well positionned) by interpolating
xyzhat_eq(1,:) = linspace(xyzhat(1,1),xyzhat(end,1),size(xyzhat,1));
xyzhat_eq(2,:) = linspace(xyzhat(1,2),xyzhat(end,2),size(xyzhat,1));
xyzhat_eq(3,:) = linspace(xyzhat(1,3),xyzhat(end,3),size(xyzhat,1));

xyzhat_eq = xyzhat_eq';

end